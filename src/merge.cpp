#include "merge.hpp"

Merge::Merge(std::string_view reference_sequence)
    : reference_sequence_{reference_sequence} {}

void Merge::AddTrees(const std::vector<std::reference_wrapper<MADAG>>& trees,
                     bool show_progress) {
  std::vector<size_t> tree_idxs;
  tree_idxs.resize(trees.size());
  std::iota(tree_idxs.begin(), tree_idxs.end(), trees_.size());

  trees_.insert(trees_.end(), trees.begin(), trees.end());
  tree_labels_.resize(trees_.size());

  ComputeCompactGenomes(tree_idxs, show_progress);
  ComputeLeafSets(tree_idxs, show_progress);
  MergeTrees(tree_idxs);
  Assert(result_nodes_.size() == result_dag_.GetNodesCount());
  Assert(result_edges_.size() == result_dag_.GetEdgesCount());
  result_dag_.BuildConnections();
}

void Merge::AddDAGs(const std::vector<std::reference_wrapper<MADAG>>& dags,
                    bool show_progress) {
  std::vector<size_t> tree_idxs;
  tree_idxs.resize(dags.size());
  std::iota(tree_idxs.begin(), tree_idxs.end(), trees_.size());

  trees_.insert(trees_.end(), dags.begin(), dags.end());
  tree_labels_.resize(trees_.size());

  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    MADAG& tree = trees_.at(tree_idx).get();
    std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    labels.resize(tree.dag.GetNodesCount());
    for (size_t node_idx = 0; node_idx < tree.dag.GetNodesCount(); ++node_idx) {
      auto cg_iter =
          all_compact_genomes_.insert(std::move(tree.compact_genomes.at(node_idx)));
      labels.at(node_idx).SetCompactGenome(std::addressof(*cg_iter.first));
    }
    tree.compact_genomes.resize(0);
  });

  ComputeLeafSets(tree_idxs, show_progress);
  MergeTrees(tree_idxs);
  Assert(result_nodes_.size() == result_dag_.GetNodesCount());
  Assert(result_edges_.size() == result_dag_.GetEdgesCount());
  result_dag_.BuildConnections();
}

DAG& Merge::GetResult() { return result_dag_; }

const DAG& Merge::GetResult() const { return result_dag_; }

const std::unordered_map<NodeLabel, NodeId>& Merge::GetResultNodes() const {
  return result_nodes_;
}

std::vector<EdgeMutations> Merge::ComputeResultEdgeMutations() const {
  std::vector<EdgeMutations> result;
  result.resize(result_dag_.GetEdgesCount());
  Assert(result_edges_.size() == result.size());
  for (auto& [label, edge_id] : result_edges_) {
    Assert(label.GetParent().GetCompactGenome());
    Assert(label.GetChild().GetCompactGenome());
    const CompactGenome& parent = *label.GetParent().GetCompactGenome();
    const CompactGenome& child = *label.GetChild().GetCompactGenome();
    EdgeMutations& muts = result.at(edge_id.value);
    muts = CompactGenome::ToEdgeMutations(reference_sequence_, parent, child);
  }
  return result;
}

void Merge::ComputeCompactGenomes(const std::vector<size_t>& tree_idxs,
                                  bool show_progress) {
  if (show_progress) {
    std::cout << "Computing compact genomes " << std::flush;
  }
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const MADAG& tree = trees_.at(tree_idx).get();
    std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    labels.resize(tree.dag.GetNodesCount());
    std::vector<CompactGenome> computed_cgs =
        tree.ComputeCompactGenomes(reference_sequence_);
    for (size_t node_idx = 0; node_idx < tree.dag.GetNodesCount(); ++node_idx) {
      auto cg_iter = all_compact_genomes_.insert(std::move(computed_cgs.at(node_idx)));
      labels.at(node_idx).SetCompactGenome(std::addressof(*cg_iter.first));
    }
    if (show_progress) {
      std::cout << "." << std::flush;
    }
  });
  if (show_progress) {
    std::cout << " done.\n";
  }
}

void Merge::ComputeLeafSets(const std::vector<size_t>& tree_idxs, bool show_progress) {
  if (show_progress) {
    std::cout << "Computing leaf sets " << std::flush;
  }
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const MADAG& tree = trees_.at(tree_idx).get();
    std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    std::vector<LeafSet> computed_ls = ComputeLeafSetsDAG(tree, labels);
    for (size_t node_idx = 0; node_idx < tree.dag.GetNodesCount(); ++node_idx) {
      auto ls_iter = all_leaf_sets_.insert(std::move(computed_ls.at(node_idx)));
      labels.at(node_idx).SetLeafSet(std::addressof(*ls_iter.first));
    }
    if (show_progress) {
      std::cout << "." << std::flush;
    }
  });
  if (show_progress) {
    std::cout << " done.\n";
  }
}

void Merge::MergeTrees(const std::vector<size_t>& tree_idxs) {
  NodeId node_id{result_dag_.GetNodesCount()};
  std::mutex mtx;
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    for (auto label : labels) {
      std::unique_lock<std::mutex> lock{mtx};
      if (result_nodes_.try_emplace(label, node_id).second) {
        ++node_id.value;
      }
    }
  });
  tbb::concurrent_vector<EdgeLabel> added_edges;
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const MADAG& tree = trees_.at(tree_idx).get();
    const std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    for (Edge edge : tree.dag.GetEdges()) {
      auto& parent_label = labels.at(edge.GetParentId().value);
      auto& child_label = labels.at(edge.GetChildId().value);
      auto ins = result_edges_.insert({{parent_label, child_label}, {}});
      if (ins.second) {
        added_edges.push_back({parent_label, child_label});
      }
    }
  });
  result_dag_.InitializeNodes(result_nodes_.size());
  EdgeId edge_id{result_dag_.GetEdgesCount()};
  for (auto edge : added_edges) {
    auto parent = result_nodes_.find(edge.GetParent());
    auto child = result_nodes_.find(edge.GetChild());
    Assert(parent != result_nodes_.end());
    Assert(child != result_nodes_.end());
    Assert(parent->second.value < result_dag_.GetNodesCount());
    Assert(child->second.value < result_dag_.GetNodesCount());
    result_dag_.AddEdge(edge_id, parent->second, child->second, {0});
    auto result_edge_it = result_edges_.find(edge);
    Assert(result_edge_it != result_edges_.end());
    result_edge_it->second = edge_id;
    edge_id.value++;
  }
}

std::vector<LeafSet> Merge::ComputeLeafSets(const MADAG& tree,
                                            const std::vector<NodeLabel>& labels) {
  std::vector<LeafSet> result;
  result.resize(tree.dag.GetNodesCount());
  for (Node node : tree.dag.TraversePostOrder()) {
    result.at(node.GetId().value) = LeafSet{node, labels, result};
  }
  return result;
}

std::vector<LeafSet> Merge::ComputeLeafSetsDAG(const MADAG& dag,
                                               const std::vector<NodeLabel>& labels) {
  std::vector<LeafSet> result;
  result.resize(dag.dag.GetNodesCount());
  auto ComputeLS = [&](auto& self, Node node) {
    LeafSet& leaf_set = result.at(node.GetId().value);
    if (not leaf_set.empty()) {
      return;
    }
    for (Node child : node.GetChildren() | Transform::GetChild()) {
      self(self, child);
    }
    result.at(node.GetId().value) = LeafSet{node, labels, result};
  };
  for (Node node : dag.dag.GetNodes()) {
    ComputeLS(ComputeLS, node);
  }
  return result;
}

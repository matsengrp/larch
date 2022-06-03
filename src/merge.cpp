#include "merge.hpp"

Merge::Merge(std::string_view reference_sequence)
    : reference_sequence_{reference_sequence} {}

void Merge::AddTrees(const std::vector<std::reference_wrapper<const MADAG>>& trees,
                     bool show_progress) {
  std::vector<size_t> tree_idxs;
  tree_idxs.resize(trees.size());
  std::iota(tree_idxs.begin(), tree_idxs.end(), trees_.size());

  trees_.insert(trees_.end(), trees.begin(), trees.end());
  tree_labels_.resize(trees_.size());

  ComputeCompactGenomes(tree_idxs, show_progress);
  ComputeLeafSets(tree_idxs, show_progress);
  MergeTrees(tree_idxs);
  Assert(result_nodes_.size() == result_dag_.GetNodes().size());
  Assert(result_edges_.size() == result_dag_.GetEdges().size());
  result_dag_.BuildConnections();
}

void Merge::AddDAGs(const std::vector<std::reference_wrapper<const MADAG>>& dags,
                    std::vector<std::vector<CompactGenome>>&& compact_genomes,
                    bool show_progress) {
  std::vector<size_t> tree_idxs;
  tree_idxs.resize(dags.size());
  std::iota(tree_idxs.begin(), tree_idxs.end(), trees_.size());

  trees_.insert(trees_.end(), dags.begin(), dags.end());
  tree_labels_.resize(trees_.size());

  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const MADAG& tree = trees_.at(tree_idx).get();
    std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    labels.resize(tree.dag.GetNodes().size());
    for (size_t node_idx = 0; node_idx < tree.dag.GetNodes().size(); ++node_idx) {
      auto cg_iter = all_compact_genomes_.insert(
          std::move(compact_genomes.at(tree_idx).at(node_idx)));
      labels.at(node_idx).compact_genome = std::addressof(*cg_iter.first);
    }
  });

  ComputeLeafSets(tree_idxs, show_progress);
  MergeTrees(tree_idxs);
  Assert(result_nodes_.size() == result_dag_.GetNodes().size());
  Assert(result_edges_.size() == result_dag_.GetEdges().size());
  result_dag_.BuildConnections();
}

DAG& Merge::GetResult() { return result_dag_; }

const DAG& Merge::GetResult() const { return result_dag_; }

const std::unordered_map<NodeLabel, NodeId>& Merge::GetResultNodes() const {
  return result_nodes_;
}

std::vector<EdgeMutations> Merge::ComputeResultEdgeMutations() const {
  std::vector<EdgeMutations> result;
  result.resize(result_dag_.GetEdges().size());
  Assert(result_edges_.size() == result.size());
  for (auto& [label, edge_id] : result_edges_) {
    Assert(label.parent_compact_genome);
    Assert(label.child_compact_genome);
    const CompactGenome& parent = *label.parent_compact_genome;
    const CompactGenome& child = *label.child_compact_genome;
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
    labels.resize(tree.dag.GetNodes().size());
    std::vector<CompactGenome> computed_cgs =
        ComputeCompactGenomesDAG(tree, reference_sequence_);
    for (size_t node_idx = 0; node_idx < tree.dag.GetNodes().size(); ++node_idx) {
      auto cg_iter = all_compact_genomes_.insert(std::move(computed_cgs.at(node_idx)));
      labels.at(node_idx).compact_genome = std::addressof(*cg_iter.first);
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
    for (size_t node_idx = 0; node_idx < tree.dag.GetNodes().size(); ++node_idx) {
      auto ls_iter = all_leaf_sets_.insert(std::move(computed_ls.at(node_idx)));
      labels.at(node_idx).leaf_set = std::addressof(*ls_iter.first);
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
  NodeId node_id{result_dag_.GetNodes().size()};
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
  tbb::concurrent_vector<std::pair<EdgeLabel, EdgeId>> added_edges;
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const MADAG& tree = trees_.at(tree_idx).get();
    const std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    EdgeId edge_id{result_dag_.GetEdges().size()};
    for (Edge edge : tree.dag.GetEdges()) {
      auto& parent_label = labels.at(edge.GetParentId().value);
      auto& child_label = labels.at(edge.GetChildId().value);
      auto ins =
          result_edges_.insert({{parent_label.compact_genome, parent_label.leaf_set,
                                 child_label.compact_genome, child_label.leaf_set},
                                {}});
      if (ins.second) {
        ins.first->second = edge_id;
        edge_id.value++;
        added_edges.push_back(*ins.first);
      }
    }
  });
  result_dag_.InitializeNodes(result_nodes_.size());
  EdgeId edge_id{result_dag_.GetEdges().size()};
  for (auto& [edge, id] : added_edges) {
    auto parent =
        result_nodes_.find(NodeLabel{edge.parent_compact_genome, edge.parent_leaf_set});
    auto child =
        result_nodes_.find(NodeLabel{edge.child_compact_genome, edge.child_leaf_set});
    Assert(parent != result_nodes_.end());
    Assert(child != result_nodes_.end());
    Assert(parent->second.value < result_dag_.GetNodes().size());
    Assert(child->second.value < result_dag_.GetNodes().size());
    id = edge_id;
    result_dag_.AddEdge(edge_id, parent->second, child->second, {0});
    edge_id.value++;
  }
}

std::vector<CompactGenome> Merge::ComputeCompactGenomes(
    const MADAG& tree, std::string_view reference_sequence) {
  std::vector<CompactGenome> result;
  result.resize(tree.dag.GetNodes().size());
  for (auto [node, edge] : tree.dag.TraversePreOrder()) {
    const EdgeMutations& mutations = tree.edge_mutations.at(edge.GetId().value);
    const CompactGenome& parent = result.at(edge.GetParentId().value);
    CompactGenome& compact_genome = result.at(node.GetId().value);
    compact_genome = CompactGenome{mutations, parent, reference_sequence};
  }
  return result;
}

std::vector<CompactGenome> Merge::ComputeCompactGenomesDAG(
    const MADAG& dag, std::string_view reference_sequence) {
  std::vector<CompactGenome> result;
  result.resize(dag.dag.GetNodes().size());
  auto ComputeCG = [&](auto& self, Node node) {
    CompactGenome& compact_genome = result.at(node.GetId().value);
    if (node.IsRoot()) {
      compact_genome = CompactGenome(node, dag.edge_mutations, reference_sequence);
      return;
    }
    if (not compact_genome.empty()) {
      return;
    }
    Edge edge = *node.GetParents().begin();
    self(self, edge.GetParent());
    const EdgeMutations& mutations = dag.edge_mutations.at(edge.GetId().value);
    const CompactGenome& parent = result.at(edge.GetParentId().value);
    compact_genome = CompactGenome{mutations, parent, reference_sequence};
  };
  for (Node node : dag.dag.GetNodes()) {
    ComputeCG(ComputeCG, node);
  }
  return result;
}

std::vector<LeafSet> Merge::ComputeLeafSets(const MADAG& tree,
                                            const std::vector<NodeLabel>& labels) {
  std::vector<LeafSet> result;
  result.resize(tree.dag.GetNodes().size());
  for (Node node : tree.dag.TraversePostOrder()) {
    result.at(node.GetId().value) = LeafSet{node, labels, result};
  }
  return result;
}

std::vector<LeafSet> Merge::ComputeLeafSetsDAG(const MADAG& dag,
                                               const std::vector<NodeLabel>& labels) {
  std::vector<LeafSet> result;
  result.resize(dag.dag.GetNodes().size());
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

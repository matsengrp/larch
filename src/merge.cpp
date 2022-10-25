#include "merge.hpp"
#include "common.hpp"

Merge::Merge(std::string_view reference_sequence) : result_dag_{reference_sequence} {}

void Merge::AddDAGs(const std::vector<std::reference_wrapper<MADAG>>& trees,
                    bool have_compact_genomes) {
  for (auto tree : trees) {
    tree.get().AssertUA();
  }
  std::vector<size_t> tree_idxs;
  tree_idxs.resize(trees.size());
  std::iota(tree_idxs.begin(), tree_idxs.end(), trees_.size());

  trees_.insert(trees_.end(), trees.begin(), trees.end());
  tree_labels_.resize(trees_.size());

  if (have_compact_genomes) {
    tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
      MADAG& tree = trees_.at(tree_idx).get();
      std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
      labels.resize(tree.GetDAG().GetNodesCount());
      for (size_t node_idx = 0; node_idx < tree.GetDAG().GetNodesCount(); ++node_idx) {
        auto cg_iter =
            all_compact_genomes_.insert(tree.ExtractCompactGenome({node_idx}));
        labels.at(node_idx).SetCompactGenome(std::addressof(*cg_iter.first));
      }
      tree.RemoveCompactGenomes();
    });
  } else {
    ComputeCompactGenomes(tree_idxs);
  }

  ComputeLeafSets(tree_idxs);
  MergeTrees(tree_idxs);
  Assert(result_nodes_.size() == result_dag_.GetDAG().GetNodesCount());
  Assert(result_edges_.size() == result_dag_.GetDAG().GetEdgesCount());
  result_dag_.BuildConnections();
  result_dag_.RemoveCompactGenomes();
  result_dag_.RemoveEdgeMutations();
}

MADAG& Merge::GetResult() { return result_dag_; }

const MADAG& Merge::GetResult() const { return result_dag_; }

const std::unordered_map<NodeLabel, NodeId>& Merge::GetResultNodes() const {
  return result_nodes_;
}

const std::vector<NodeLabel>& Merge::GetResultNodeLabels() const {
  return result_node_labels_;
}

void Merge::ComputeResultEdgeMutations() {
  std::vector<EdgeMutations> result;
  result.resize(result_dag_.GetDAG().GetEdgesCount());
  Assert(result_edges_.size() == result.size());
  for (auto& [label, edge_id] : result_edges_) {
    Assert(label.GetParent().GetCompactGenome());
    Assert(label.GetChild().GetCompactGenome());
    const CompactGenome& parent = *label.GetParent().GetCompactGenome();
    const CompactGenome& child = *label.GetChild().GetCompactGenome();
    EdgeMutations& muts = result.at(edge_id.value);
    muts = CompactGenome::ToEdgeMutations(result_dag_.GetReferenceSequence(), parent,
                                          child);
  }
  result_dag_.SetEdgeMutations(std::move(result));
}

bool Merge::ContainsLeafset(const LeafSet& leafset) const {
  return all_leaf_sets_.find(leafset) != all_leaf_sets_.end();
}

void Merge::ComputeCompactGenomes(const std::vector<size_t>& tree_idxs) {
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const MADAG& tree = trees_.at(tree_idx).get();
    std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    labels.resize(tree.GetDAG().GetNodesCount());
    std::vector<CompactGenome> computed_cgs = tree.ComputeCompactGenomes();
    for (size_t node_idx = 0; node_idx < tree.GetDAG().GetNodesCount(); ++node_idx) {
      auto cg_iter = all_compact_genomes_.insert(std::move(computed_cgs.at(node_idx)));
      labels.at(node_idx).SetCompactGenome(std::addressof(*cg_iter.first));
    }
  });
}

void Merge::ComputeLeafSets(const std::vector<size_t>& tree_idxs) {
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const MADAG& tree = trees_.at(tree_idx).get();
    std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    std::vector<LeafSet> computed_ls = ComputeLeafSets(tree, labels);
    for (size_t node_idx = 0; node_idx < tree.GetDAG().GetNodesCount(); ++node_idx) {
      auto ls_iter = all_leaf_sets_.insert(std::move(computed_ls.at(node_idx)));
      labels.at(node_idx).SetLeafSet(std::addressof(*ls_iter.first));
    }
  });
}

void Merge::MergeTrees(const std::vector<size_t>& tree_idxs) {
  NodeId node_id{result_dag_.GetDAG().GetNodesCount()};
  std::mutex mtx;
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    for (auto label : labels) {
      std::unique_lock<std::mutex> lock{mtx};
      if (result_nodes_.try_emplace(label, node_id).second) {
        GetOrInsert(result_node_labels_, node_id) = label;
        ++node_id.value;
      }
    }
  });
  tbb::concurrent_vector<EdgeLabel> added_edges;
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const MADAG& tree = trees_.at(tree_idx).get();
    const std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    for (Edge edge : tree.GetDAG().GetEdges()) {
      auto& parent_label = labels.at(edge.GetParentId().value);
      auto& child_label = labels.at(edge.GetChildId().value);
      auto ins = result_edges_.insert({{parent_label, child_label}, {}});
      if (ins.second) {
        added_edges.push_back({parent_label, child_label});
      }
    }
  });
  result_dag_.InitializeNodes(result_nodes_.size());
  EdgeId edge_id{result_dag_.GetDAG().GetEdgesCount()};
  for (auto edge : added_edges) {
    auto parent = result_nodes_.find(edge.GetParent());
    auto child = result_nodes_.find(edge.GetChild());
    Assert(parent != result_nodes_.end());
    Assert(child != result_nodes_.end());
    Assert(parent->second.value < result_dag_.GetDAG().GetNodesCount());
    Assert(child->second.value < result_dag_.GetDAG().GetNodesCount());
    result_dag_.AddEdge(edge_id, parent->second, child->second, edge.ComputeCladeIdx());
    auto result_edge_it = result_edges_.find(edge);
    Assert(result_edge_it != result_edges_.end());
    result_edge_it->second = edge_id;
    edge_id.value++;
  }
}

std::vector<LeafSet> Merge::ComputeLeafSets(const MADAG& dag,
                                            const std::vector<NodeLabel>& labels) {
  std::vector<LeafSet> result;
  result.resize(dag.GetDAG().GetNodesCount());
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
  for (Node node : dag.GetDAG().GetNodes()) {
    ComputeLS(ComputeLS, node);
  }
  return result;
}

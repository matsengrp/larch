#include <mutex>

Merge::Merge(std::string_view reference_sequence) {
  ResultDAG().SetReferenceSequence(reference_sequence);
}

template <typename DAGSRange>
void Merge::AddDAGs(const DAGSRange& dags, NodeId below) {
  if (below.value != NoId and dags.size() != 1) {
    Fail("Pass exactly one DAG when merging below given node.");
  }
  if (below.value != NoId and dags.at(0).Get(below).IsUA()) {
    below.value = NoId;
  }

  std::vector<size_t> idxs;
  idxs.resize(dags.size());
  std::iota(idxs.begin(), idxs.end(), 0);

  std::vector<std::vector<NodeLabel>> dags_labels;
  dags_labels.resize(dags.size());

  tbb::parallel_for_each(idxs, [&](size_t i) {
    auto dag = dags.at(i).GetRoot().GetDAG();
    dag.AssertUA();
    auto& labels = dags_labels.at(i);
    labels.resize(dag.GetNodesCount());
    for (auto node : dag.GetNodes()) {
      if (below.value != NoId and node.IsUA()) {
        continue;
      }
      auto cg_iter =
          ResultDAG().AddDeduplicated(node.Const().GetCompactGenome().Copy());
      labels.at(node.GetId().value).SetCompactGenome(cg_iter.first);
    }
  });

  ComputeLeafSets(dags, dags_labels, below);
  MergeDAGs(dags, dags_labels, below);
  Assert(result_nodes_.size() == ResultDAG().GetNodesCount());
  Assert(result_node_labels_.size() == ResultDAG().GetNodesCount());
  Assert(result_edges_.size() == ResultDAG().GetEdgesCount());
  ResultDAG().BuildConnections();
  ComputeResultEdgeMutations();
}

MergeDAG Merge::GetResult() const { return result_dag_storage_.View(); }

MutableMergeDAG Merge::ResultDAG() { return result_dag_storage_.View(); }

const std::unordered_map<NodeLabel, NodeId>& Merge::GetResultNodes() const {
  return result_nodes_;
}

const std::vector<NodeLabel>& Merge::GetResultNodeLabels() const {
  return result_node_labels_;
}

void Merge::ComputeResultEdgeMutations() {
  for (auto& [label, edge_id] : result_edges_) {
    Assert(label.GetParent().GetCompactGenome());
    Assert(label.GetChild().GetCompactGenome());
    const CompactGenome& parent = *label.GetParent().GetCompactGenome();
    const CompactGenome& child = *label.GetChild().GetCompactGenome();
    ResultDAG().Get(edge_id).SetEdgeMutations(CompactGenome::ToEdgeMutations(
        ResultDAG().GetReferenceSequence(), parent, child));
  }
}

bool Merge::ContainsLeafset(const LeafSet& leafset) const {
  return all_leaf_sets_.find(leafset) != all_leaf_sets_.end();
}

template <typename DAGSRange>
void Merge::ComputeLeafSets(const DAGSRange& dags,
                            std::vector<std::vector<NodeLabel>>& dags_labels,
                            NodeId below) {
  std::vector<size_t> idxs;
  idxs.resize(dags.size());
  std::iota(idxs.begin(), idxs.end(), 0);
  tbb::parallel_for_each(idxs, [&](size_t i) {
    auto dag = dags.at(i).GetRoot().GetDAG();
    std::vector<NodeLabel>& labels = dags_labels.at(i);
    std::vector<LeafSet> computed_ls = LeafSet::ComputeLeafSets(dag, labels);
    for (auto node : dag.GetNodes()) {
      if (below.value != NoId and node.IsUA()) {
        continue;
      }
      auto ls_iter =
          all_leaf_sets_.insert(std::move(computed_ls.at(node.GetId().value)));
      labels.at(node.GetId().value).SetLeafSet(std::addressof(*ls_iter.first));
    }
  });
}

template <typename DAGSRange>
void Merge::MergeDAGs(const DAGSRange& dags,
                      const std::vector<std::vector<NodeLabel>>& dags_labels,
                      NodeId below) {
  NodeId node_id{ResultDAG().GetNodesCount()};
  std::vector<size_t> idxs;
  idxs.resize(dags.size());
  std::iota(idxs.begin(), idxs.end(), 0);
  std::mutex mtx;
  tbb::parallel_for_each(idxs, [&](size_t idx) {
    NodeId id{0};
    for (auto& label : dags_labels.at(idx)) {
      label.AssertNonEmpty();
      auto dag = dags.at(idx);
      std::unique_lock<std::mutex> lock{mtx};
      auto insert_pair = result_nodes_.insert({label, node_id});
      if (insert_pair.second) {
        GetOrInsert(result_node_labels_, node_id) = label;
        if constexpr (decltype(dag)::template contains_element_feature<NodeId,
                                                                       MappedNodes>) {
          dag.Get(id).SetOriginalId(node_id);
        }
        ++node_id.value;
      } else {
        if (id.value != dag.GetNodesCount() - 1) {
          if constexpr (decltype(dag)::template contains_element_feature<NodeId,
                                                                         MappedNodes>) {
            dag.Get(id).SetOriginalId(insert_pair.first->second);
          }
        }
      }
      ++id.value;
    }
  });
  tbb::concurrent_vector<EdgeLabel> added_edges;
  tbb::parallel_for_each(idxs, [&](size_t idx) {
    auto& dag = dags.at(idx);
    const std::vector<NodeLabel>& labels = dags_labels.at(idx);
    for (auto edge : dag.GetEdges()) {
      if (below.value != NoId and edge.IsUA()) {
        continue;
      }
      const auto& parent_label = labels.at(edge.GetParentId().value);
      parent_label.AssertNonEmpty();
      const auto& child_label = labels.at(edge.GetChildId().value);
      child_label.AssertNonEmpty();
      auto ins = result_edges_.insert({{parent_label, child_label}, {}});
      if (ins.second) {
        added_edges.push_back({parent_label, child_label});
      }
    }
  });
  if (below.value != NoId) {
    auto below_node = dags.at(0).Get(below);
    EdgeLabel below_edge = {
        result_node_labels_.at(below_node.GetFirstParent().GetParent().GetId().value),
        dags_labels.at(0).at(
            dags.at(0).GetRoot().GetFirstChild().GetChild().GetId().value)};
    auto ins = result_edges_.insert({below_edge, {}});
    if (ins.second) {
      added_edges.push_back(below_edge);
    }
  }
  ResultDAG().InitializeNodes(result_nodes_.size());
  EdgeId edge_id{ResultDAG().GetEdgesCount()};
  for (auto& edge : added_edges) {
    auto parent = result_nodes_.find(edge.GetParent());
    auto child = result_nodes_.find(edge.GetChild());
    Assert(parent != result_nodes_.end());
    Assert(child != result_nodes_.end());
    Assert(parent->second.value < ResultDAG().GetNodesCount());
    Assert(child->second.value < ResultDAG().GetNodesCount());
    CladeIdx clade = edge.ComputeCladeIdx();
    Assert(clade.value != NoId);
    ResultDAG().AddEdge(edge_id, parent->second, child->second, clade);
    ResultDAG().Get(parent->second) = parent->first.GetCompactGenome();
    ResultDAG().Get(child->second) = child->first.GetCompactGenome();
    auto result_edge_it = result_edges_.find(edge);
    Assert(result_edge_it != result_edges_.end());
    result_edge_it->second = edge_id;
    edge_id.value++;
  }
}

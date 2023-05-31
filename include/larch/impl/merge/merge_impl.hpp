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
    for (auto node : dag.Const().GetNodes()) {
      if (below.value != NoId and node.IsUA()) {
        continue;
      }
      auto cg_iter = ResultDAG().AddDeduplicated(node.GetCompactGenome().Copy());
      labels.at(node.GetId().value).SetCompactGenome(cg_iter.first);
    }
  });

  tbb::parallel_for_each(idxs, [&](size_t i) {
    auto dag = dags.at(i).GetRoot().GetDAG();
    std::vector<NodeLabel>& labels = dags_labels.at(i);
    std::vector<LeafSet> computed_ls = LeafSet::ComputeLeafSets(dag, labels);
    for (auto node : dag.GetNodes()) {
      if (below.value != NoId and node.IsUA()) {
        continue;
      }
      auto& ls = computed_ls.at(node.GetId().value);
      auto& label = labels.at(node.GetId().value);
      auto ls_iter = all_leaf_sets_.insert(std::move(ls));
      label.SetLeafSet(std::addressof(*ls_iter.first));
      Assert(not label.Empty());
    }
  });

  NodeId node_id{ResultDAG().GetNodesCount()};
  std::mutex mtx;
  tbb::parallel_for_each(idxs, [&](size_t idx) {
    NodeId id{0};
    auto& dag = dags.at(idx);
    auto& labels = dags_labels.at(idx);
    for (auto node : dag.GetNodes()) {
      auto& label = labels.at(node.GetId().value);
      if (below.value != NoId and node.IsUA()) {
        continue;
      }
      Assert(not label.Empty());
      auto [insert_pair, orig_id] = [&] {
        std::unique_lock<std::mutex> lock{mtx};
        auto ins_pair = result_nodes_.insert({label, node_id});
        if (ins_pair.second) {
          GetOrInsert(result_node_labels_, node_id) = label;
        }
        auto result = std::make_pair(ins_pair, node_id);
        if (ins_pair.second) {
          ++node_id.value;
        }
        return result;
      }();
      if (insert_pair.second) {
        if constexpr (std::decay_t<decltype(dag)>::template contains_element_feature<
                          NodeId, MappedNodes>) {
          dag.Get(id).SetOriginalId(orig_id);
        }
      } else {
        if (id.value != dag.GetNodesCount() - 1) {
          if constexpr (std::decay_t<decltype(dag)>::template contains_element_feature<
                            NodeId, MappedNodes>) {
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
      const auto& child_label = labels.at(edge.GetChildId().value);
      Assert(not parent_label.Empty());
      Assert(not child_label.Empty());
      auto ins = result_edges_.insert({{parent_label, child_label}, {}});
      if (ins.second) {
        added_edges.push_back({parent_label, child_label});
      }
    }
  });

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

  Assert(result_nodes_.size() == ResultDAG().GetNodesCount());
  Assert(result_node_labels_.size() == ResultDAG().GetNodesCount());
  Assert(result_edges_.size() == ResultDAG().GetEdgesCount());
  ResultDAG().BuildConnections();
  ComputeResultEdgeMutations();
}

template <typename DAG>
void Merge::AddDAG(DAG& dag, NodeId below) {
  if (below.value != NoId and dag.Get(below).IsUA()) {
    below.value = NoId;
  }
  struct {
    DAG& dag_;
    size_t size() const { return 1; }
    const auto* begin() const { return &dag_; }
    const auto* end() const { return &dag_ + 1; }
    const auto& at(size_t i) const {
      Assert(i == 0);
      return dag_;
    }
    auto* begin() { return &dag_; }
    auto* end() { return &dag_ + 1; }
    auto& at(size_t i) {
      Assert(i == 0);
      return dag_;
    }
  } dags{dag};
  AddDAGs(dags, below);
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

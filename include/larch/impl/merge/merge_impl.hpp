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

  if (dags.size() == 0) {
    return;
  }

  const bool was_empty = ResultDAG().IsEmpty();

  std::vector<size_t> idxs;
  idxs.resize(dags.size());
  std::iota(idxs.begin(), idxs.end(), 0);

  std::vector<std::vector<NodeLabel>> dags_labels;
  dags_labels.resize(dags.size());

  tbb::parallel_for_each(idxs, [&](size_t i) {
    MergeCompactGenomes(i, dags, below, dags_labels, ResultDAG());
  });

  tbb::parallel_for_each(idxs, [&](size_t i) {
    ComputeLeafSets(i, dags, below, dags_labels, all_leaf_sets_);
  });

  std::atomic<size_t> node_id{ResultDAG().GetNodesCount()};
  tbb::concurrent_vector<std::tuple<EdgeLabel, EdgeId, NodeId, NodeId, CladeIdx>>
      added_edges;
  tbb::parallel_for_each(idxs, [&](size_t i) {
    MergeNodes(i, dags, below, dags_labels, result_nodes_, result_node_labels_,
               node_id);
    MergeEdges(i, dags, below, dags_labels, result_edges_, added_edges);
  });

  ResultDAG().InitializeNodes(result_nodes_.size());
  std::atomic<size_t> edge_id{ResultDAG().GetEdgesCount()};
  ResultDAG().InitializeEdges(result_edges_.size());
  idxs.resize(added_edges.size());
  std::iota(idxs.begin(), idxs.end(), 0);
  tbb::parallel_for_each(idxs, [&](size_t i) {
    BuildResult(i, added_edges, edge_id, result_nodes_, result_edges_, ResultDAG());
  });

  if (was_empty) {
    ResultDAG().BuildConnections();
  } else {
    for (auto& [label, id, parent_id, child_id, clade] : added_edges) {
      ResultDAG().Get(parent_id).AddEdge(clade, id, true);
      ResultDAG().Get(child_id).AddEdge(clade, id, false);
    }
  }

  Assert(result_nodes_.size() == ResultDAG().GetNodesCount());
  Assert(result_node_labels_.size() == ResultDAG().GetNodesCount());
  Assert(result_edges_.size() == ResultDAG().GetEdgesCount());
}

template <typename DAG>
void Merge::AddDAG(DAG& dag, NodeId below) {
  if (below.value != NoId and dag.Get(below).IsUA()) {
    below.value = NoId;
  }
  struct {
    DAG& dag_;  // NOLINT
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

const ConcurrentUnorderedMap<NodeLabel, NodeId>& Merge::GetResultNodes() const {
  return result_nodes_;
}

const ConcurrentUnorderedMap<NodeId, NodeLabel>& Merge::GetResultNodeLabels() const {
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
void Merge::MergeCompactGenomes(size_t i, const DAGSRange& dags, NodeId below,
                                std::vector<std::vector<NodeLabel>>& dags_labels,
                                MutableMergeDAG result_dag) {
  auto dag = dags.at(i).GetRoot().GetDAG();
  dag.AssertUA();
  auto& labels = dags_labels.at(i);
  labels.resize(dag.GetNodesCount());
  for (auto node : dag.Const().GetNodes()) {
    if (below.value != NoId and node.IsUA()) {
      continue;
    }
    auto cg_iter = result_dag.AddDeduplicated(node.GetCompactGenome());
    labels.at(node.GetId().value).SetCompactGenome(cg_iter.first);
  }
}

template <typename DAGSRange>
void Merge::ComputeLeafSets(size_t i, const DAGSRange& dags, NodeId below,
                            std::vector<std::vector<NodeLabel>>& dags_labels,
                            ConcurrentUnorderedSet<LeafSet>& all_leaf_sets) {
  auto dag = dags.at(i).GetRoot().GetDAG();
  std::vector<NodeLabel>& labels = dags_labels.at(i);
  std::vector<LeafSet> computed_ls = LeafSet::ComputeLeafSets(dag, labels);
  for (auto node : dag.GetNodes()) {
    if (below.value != NoId and node.IsUA()) {
      continue;
    }
    auto& ls = computed_ls.at(node.GetId().value);
    auto& label = labels.at(node.GetId().value);
    auto ls_iter = all_leaf_sets.insert(std::move(ls));
    label.SetLeafSet(std::addressof(*ls_iter.first));
    Assert(not label.Empty());
  }
}

template <typename DAGSRange>
void Merge::MergeNodes(size_t i, const DAGSRange& dags, NodeId below,
                       std::vector<std::vector<NodeLabel>>& dags_labels,
                       ConcurrentUnorderedMap<NodeLabel, NodeId>& result_nodes,
                       ConcurrentUnorderedMap<NodeId, NodeLabel>& result_node_labels,
                       std::atomic<size_t>& node_id) {
  NodeId id{0};
  auto& dag = dags.at(i);
  auto& labels = dags_labels.at(i);
  for (auto node : dag.GetNodes()) {
    auto& label = labels.at(node.GetId().value);
    if (below.value != NoId and node.IsUA()) {
      continue;
    }
    Assert(not label.Empty());
    auto [insert_pair, orig_id] = [&] {
      NodeId new_id;
      auto ins_pair = result_nodes.insert({label, new_id});
      if (ins_pair.second) {
        new_id.value = node_id.fetch_add(1);
        ins_pair.first->second = new_id;
        result_node_labels[new_id] = label;
      } else {
        new_id.value = ins_pair.first->second.value;
      }
      auto result = std::make_pair(ins_pair, new_id);
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
}

template <typename DAGSRange>
void Merge::MergeEdges(
    size_t i, const DAGSRange& dags, NodeId below,
    std::vector<std::vector<NodeLabel>>& dags_labels,
    ConcurrentUnorderedMap<EdgeLabel, EdgeId>& result_edges,
    tbb::concurrent_vector<std::tuple<EdgeLabel, EdgeId, NodeId, NodeId, CladeIdx>>&
        added_edges) {
  auto& dag = dags.at(i);
  const std::vector<NodeLabel>& labels = dags_labels.at(i);
  for (auto edge : dag.GetEdges()) {
    if (below.value != NoId and edge.IsUA()) {
      continue;
    }
    const auto& parent_label = labels.at(edge.GetParentId().value);
    const auto& child_label = labels.at(edge.GetChildId().value);
    Assert(not parent_label.Empty());
    Assert(not child_label.Empty());
    auto ins = result_edges.insert({{parent_label, child_label}, {}});
    if (ins.second) {
      added_edges.push_back({{parent_label, child_label}, {}, {}, {}, {}});
    }
  }
}

void Merge::BuildResult(
    size_t i,
    tbb::concurrent_vector<std::tuple<EdgeLabel, EdgeId, NodeId, NodeId, CladeIdx>>&
        added_edges,
    std::atomic<size_t>& edge_id,
    const ConcurrentUnorderedMap<NodeLabel, NodeId>& result_nodes,
    const ConcurrentUnorderedMap<EdgeLabel, EdgeId>& result_edges,
    MutableMergeDAG result_dag) {
  auto& [edge, id, parent_id, child_id, clade] = added_edges.at(i);
  id = {edge_id.fetch_add(1)};
  auto parent = result_nodes.find(edge.GetParent());
  auto child = result_nodes.find(edge.GetChild());
  Assert(parent != result_nodes.end());
  Assert(child != result_nodes.end());
  Assert(parent->second.value < result_dag.GetNodesCount());
  Assert(child->second.value < result_dag.GetNodesCount());
  parent_id = parent->second;
  child_id = child->second;
  clade = edge.ComputeCladeIdx();
  Assert(clade.value != NoId);
  result_dag.Get(id).Set(parent->second, child->second, clade);
  result_dag.Get(parent->second) = parent->first.GetCompactGenome();
  result_dag.Get(child->second) = child->first.GetCompactGenome();
  auto result_edge_it = result_edges.find(edge);
  Assert(result_edge_it != result_edges.end());
  result_edge_it->second = id;
}

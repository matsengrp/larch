Merge::Merge(std::string_view reference_sequence) : result_dag_storage_{{}} {
  ResultDAG().SetReferenceSequence(reference_sequence);
}

template <typename DAGSRange>
void Merge::AddDAGs(DAGSRange&& dags, NodeId below) {
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

  std::vector<std::vector<NodeLabel>> dags_labels;
  dags_labels.resize(dags.size());

  parallel_for_each(dags.size(), [&](size_t i, size_t) {
    MergeCompactGenomes(dags.at(i), below, dags_labels.at(i));
  });

  parallel_for_each(dags.size(), [&](size_t i, size_t) {
    ComputeLeafSets(dags.at(i), below, dags_labels.at(i));
  });

  std::atomic<size_t> node_id{ResultDAG().GetNodesCount()};
  Reduction<NodeId> added_nodes{DefaultScheduler().WorkersCount()};
  parallel_for_each(dags.size(), [&](size_t i, size_t worker) {
    MergeNodes(worker, dags.at(i), below, dags_labels.at(i), node_id, added_nodes);
  });

  Reduction<std::tuple<EdgeLabel, EdgeId, NodeId, NodeId, CladeIdx>> added_edges{
      DefaultScheduler().WorkersCount()};
  parallel_for_each(dags.size(), [&](size_t i, size_t worker) {
    MergeEdges(worker, dags.at(i), below, dags_labels.at(i), added_edges);
  });

  ResultDAG().InitializeNodes(result_nodes_.Size());
  std::atomic<size_t> edge_id{ResultDAG().GetEdgesCount()};
  ResultDAG().InitializeEdges(result_edges_.Size());

  for (auto& batch : added_edges.Get()) {
    parallel_for_each(batch.size(),
                      [&](size_t i, size_t) { BuildResult(batch.at(i), edge_id); });
  }

  for (auto& batch : added_nodes.Get()) {
    parallel_for_each(batch.size(), [&](size_t i, size_t) {
      NodeId id = batch.at(i);
      result_node_labels_.At(id).Get(
          [&](auto& val) { ResultDAG().Get(id) = val.GetCompactGenome(); });
    });
  }

  if (was_empty) {
    ResultDAG().BuildConnections();
  } else {
    for (const auto& [label, id, parent_id, child_id, clade] : added_edges.GetAll()) {
      ResultDAG().Get(parent_id).AddEdge(clade, id, true);
      ResultDAG().Get(child_id).AddEdge(clade, id, false);
    }
  }

  Assert(result_nodes_.Size() == ResultDAG().GetNodesCount());
  Assert(result_node_labels_.Size() == ResultDAG().GetNodesCount());
  Assert(result_edges_.Size() == ResultDAG().GetEdgesCount());
  GetResult().GetRoot().Validate(true, true);
}

template <typename DAG>
void Merge::AddDAG(DAG&& dag, NodeId below) {
  if (below.value != NoId and dag.Get(below).IsUA()) {
    below.value = NoId;
  }
  AddDAGs(ranges::views::single(dag), below);
}

MergeDAG Merge::GetResult() const { return result_dag_storage_.View(); }

const ConcurrentUnorderedMap<NodeLabel, NodeId>& Merge::GetResultNodes() const {
  return result_nodes_;
}

const ConcurrentUnorderedMap<NodeId, NodeLabel>& Merge::GetResultNodeLabels() const {
  return result_node_labels_;
}

void Merge::ComputeResultEdgeMutations() {
  for (auto& [label, edge_id] : result_edges_.All()) {
    Assert(label.GetParent().GetCompactGenome());
    Assert(label.GetChild().GetCompactGenome());
    const CompactGenome& parent = *label.GetParent().GetCompactGenome();
    const CompactGenome& child = *label.GetChild().GetCompactGenome();
    ResultDAG().Get(edge_id).SetEdgeMutations(CompactGenome::ToEdgeMutations(
        ResultDAG().GetReferenceSequence(), parent, child));
  }
}

bool Merge::ContainsLeafset(const LeafSet& leafset) const {
  return all_leaf_sets_.Contains(leafset);
}

MutableMergeDAG Merge::ResultDAG() { return result_dag_storage_.View(); }

template <typename DAG>
void Merge::MergeCompactGenomes(DAG dag, NodeId below, std::vector<NodeLabel>& labels) {
  auto full_dag = dag.GetRoot().GetDAG();
  full_dag.AssertUA();
  labels.resize(full_dag.GetNodesCount());
  for (auto node : full_dag.Const().GetNodes()) {
    if (below.value != NoId and node.IsUA()) {
      continue;
    }
    auto cg_iter = ResultDAG().AddDeduplicated(node.GetCompactGenome());
    labels.at(node.GetId().value).SetCompactGenome(cg_iter.first);
  }
}

template <typename DAG>
void Merge::ComputeLeafSets(DAG dag, NodeId below, std::vector<NodeLabel>& labels) {
  auto full_dag = dag.GetRoot().GetDAG();
  std::vector<LeafSet> computed_ls = LeafSet::ComputeLeafSets(full_dag, labels);
  for (auto node : full_dag.GetNodes()) {
    if (below.value != NoId and node.IsUA()) {
      continue;
    }
    auto& ls = computed_ls.at(node.GetId().value);
    auto& label = labels.at(node.GetId().value);
    auto ls_iter = all_leaf_sets_.Insert(std::move(ls));
    ls_iter.first.Get([&](auto& val) { label.SetLeafSet(std::addressof(val)); });
    Assert(not label.Empty());
  }
}

template <typename DAG>
void Merge::MergeNodes(size_t worker, DAG dag, NodeId below,
                       const std::vector<NodeLabel>& labels,
                       std::atomic<size_t>& node_id, Reduction<NodeId>& added_nodes) {
  NodeId id{0};
  for (auto node : dag.GetNodes()) {
    auto& label = labels.at(node.GetId().value);
    if (below.value != NoId and node.IsUA()) {
      continue;
    }
    Assert(not label.Empty());
    auto [insert_pair, orig_id] = [&] {
      NodeId new_id;
      auto ins_pair = result_nodes_.Insert({label, new_id});
      if (ins_pair.second) {
        new_id.value = node_id.fetch_add(1);
        ins_pair.first.Get([&](auto& val) { val = new_id; });
        result_node_labels_.Insert({new_id, label});
        added_nodes.Emplace(worker, new_id);
      } else {
        ins_pair.first.Get([&](auto& val) { new_id.value = val.value; });
      }
      auto result = std::make_pair(ins_pair, new_id);
      return result;
    }();
    if (insert_pair.second) {
      if constexpr (std::decay_t<decltype(dag)>::template contains_element_feature<
                        Component::Node, MappedNodes>) {
        dag.Get(id).SetOriginalId(orig_id);
      }
    } else {
      if (id.value != dag.GetNodesCount() - 1) {
        if constexpr (std::decay_t<decltype(dag)>::template contains_element_feature<
                          Component::Node, MappedNodes>) {
          insert_pair.first.Get([&](auto& val) { dag.Get(id).SetOriginalId(val); });
        }
      }
    }
    ++id.value;
  }
}

template <typename DAG>
void Merge::MergeEdges(
    size_t worker, DAG dag, NodeId below, const std::vector<NodeLabel>& labels,
    Reduction<std::tuple<EdgeLabel, EdgeId, NodeId, NodeId, CladeIdx>>& added_edges) {
  for (auto edge : dag.GetEdges()) {
    if (below.value != NoId and edge.IsUA()) {
      continue;
    }
    const auto& parent_label = labels.at(edge.GetParentId().value);
    const auto& child_label = labels.at(edge.GetChildId().value);
    Assert(not parent_label.Empty());
    Assert(not child_label.Empty());
    auto ins = result_edges_.Insert({{parent_label, child_label}, {}});
    if (ins.second) {
      added_edges.Emplace(worker, EdgeLabel{parent_label, child_label}, EdgeId{},
                          NodeId{}, NodeId{}, CladeIdx{});
    }
  }
}

void Merge::BuildResult(AddedEdge& added_edge, std::atomic<size_t>& edge_id) {
  auto& [edge, id, parent_id, child_id, clade] = added_edge;
  id = {edge_id.fetch_add(1)};
  auto parent = result_nodes_.Find(edge.GetParent());
  auto child = result_nodes_.Find(edge.GetChild());
  Assert(parent.has_value());
  Assert(child.has_value());
  parent.value().Get([&](auto& val) { parent_id = val; });
  child.value().Get([&](auto& val) { child_id = val; });
  Assert(parent_id.value < ResultDAG().GetNodesCount());
  Assert(child_id.value < ResultDAG().GetNodesCount());
  clade = edge.ComputeCladeIdx();
  Assert(clade.value != NoId);
  ResultDAG().Get(id).Set(parent_id, child_id, clade);
  auto result_edge_it = result_edges_.Find(edge);
  Assert(result_edge_it.has_value());
  result_edge_it.value().Get([&](auto& val) { val = id; });
}

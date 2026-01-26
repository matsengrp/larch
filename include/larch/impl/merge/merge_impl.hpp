
Merge::Merge(std::string_view reference_sequence)
    : result_dag_storage_{MergeDAGStorage<>::EmptyDefault()} {
  ResultDAG().SetReferenceSequence(reference_sequence);
}

template <typename DAGSRange>
void Merge::AddDAGs(const DAGSRange& dags, NodeId below) {
  std::unique_lock lock{add_dags_mtx_};
  if (below.value != NoId and dags.size() != 1) {
    Fail("Pass exactly one DAG when merging below given node.");
  }
  if (below.value != NoId and dags.at(0).Get(below).IsUA()) {
    below.value = NoId;
  }

  if (dags.size() == 0) {
    return;
  }

  const bool was_empty = ResultDAG().empty();

  std::vector<size_t> idxs;
  idxs.resize(dags.size());
  std::iota(idxs.begin(), idxs.end(), 0);

#ifdef KEEP_ASSERTS
  // TODO uncomment once extra edges/nodes are cleared in CollapseEmptyFragmentEdges
  // ParallelForEach(idxs, [&](size_t i) {
  //   auto dag = GetFullDAG(dags.at(i));
  //   dag.GetRoot().Validate(true, not dag.IsTree());
  // });
#endif

  constexpr IdContinuity id_continuity = IdContinuity::Sparse;
  // FIXME constexpr IdContinuity id_continuity =
  // std::remove_reference_t<decltype(dags.at(
  //     0))>::template id_continuity<Component::Node>;
  std::vector<IdContainer<NodeId, NodeLabel, id_continuity, Ordering::Ordered>>
      dags_labels;
  dags_labels.resize(dags.size());

  ParallelForEach(idxs,
                  [&](size_t i) { MergeCompactGenomes(i, dags, below, dags_labels); });

  SeqForEach(idxs, [&](size_t i) {
    ComputeLeafSets(i, dags, below, dags_labels);
  });  // FIXME Parallel

#ifdef KEEP_ASSERTS
  ParallelForEach(idxs, [&](size_t i) {
    for (auto node : dags.at(i).GetNodes()) {
      Assert(not dags_labels.at(i).at(node).empty());
    }
  });
#endif

  std::atomic<size_t> node_id{ResultDAG().GetNextAvailableNodeId<MergeDAG>().value};
  ParallelForEach(idxs,
                  [&](size_t i) { MergeNodes(i, dags, below, dags_labels, node_id); });

#ifdef KEEP_ASSERTS
  result_nodes_.ReadAll([](auto result_nodes) {
    for (auto& i : result_nodes) {
      Assert(i.second.value != NoId);
    }
  });
#endif

  Reduction<std::vector<AddedEdge>> added_edges_reduction{32};
  ParallelForEach(idxs, [&](size_t i) {
    MergeEdges(i, dags, below, dags_labels, added_edges_reduction);
  });
  std::vector<AddedEdge> added_edges;
  added_edges_reduction.GatherAndClear(
      [&added_edges_reduction](auto buckets, auto& result) {
        for (auto&& bucket : buckets) {
          result.reserve(added_edges_reduction.size_approx());
          result.insert(result.end(), bucket.begin(), bucket.end());
        }
      },
      added_edges);

  ResultDAG().InitializeNodes(result_nodes_.size());
  std::atomic<size_t> edge_id{
      ResultDAG().template GetNextAvailableEdgeId<MergeDAG>().value};
  ResultDAG().InitializeEdges(result_edges_.size());
  idxs.resize(added_edges.size());
  std::iota(idxs.begin(), idxs.end(), 0);
  SeqForEach(  // FIXME ParallelForEach
      idxs, [&](size_t i) { BuildResult(i, added_edges, edge_id); });  // FIXME parallel

  if (was_empty) {
    ResultDAG().BuildConnections();
  } else {
    for (auto& [label, id, parent_id, child_id, clade] : added_edges) {
      ResultDAG().Get(parent_id).AddEdge(clade, id, true);
      ResultDAG().Get(child_id).AddEdge(clade, id, false);
    }
    for ([[maybe_unused]] auto& [label, id, parent_id, child_id, clade] : added_edges) {
      if (ResultDAG().Get(child_id).IsLeaf()) {
        ResultDAG().AddLeaf(child_id);
      }
    }
  }
  Assert(result_nodes_.size() == ResultDAG().GetNodesCount());
  Assert(result_node_labels_.size() == ResultDAG().GetNodesCount());
  Assert(result_edges_.size() == ResultDAG().GetEdgesCount());
  GetResult().GetRoot().Validate(true, true);
}

template <typename DAG>
void Merge::AddDAG(DAG dag, NodeId below) {
  if (below.value != NoId and dag.Get(below).IsUA()) {
    below.value = NoId;
  }
  struct {
    DAG& dag_;  // NOLINT
    size_t size() const { return 1; }
    const auto* begin() const { return &dag_; }
    const auto* end() const { return &dag_ + 1; }
    const auto& at([[maybe_unused]] size_t i) const {
      Assert(i == 0);
      return dag_;
    }
    auto* begin() { return &dag_; }
    auto* end() { return &dag_ + 1; }
    auto& at([[maybe_unused]] size_t i) {
      Assert(i == 0);
      return dag_;
    }
  } dags{dag};
  AddDAGs(dags, below);
}

MergeDAG Merge::GetResult() const { return result_dag_storage_.View().Const(); }

MutableMergeDAG Merge::ResultDAG() { return result_dag_storage_.View(); }

const GrowableHashMap<NodeLabel, NodeId>& Merge::GetResultNodes() const {
  return result_nodes_;
}

const GrowableHashMap<NodeId, NodeLabel>& Merge::GetResultNodeLabels() const {
  return result_node_labels_;
}

const GrowableHashMap<std::string, CompactGenome>& Merge::SampleIdToCGMap() const {
  return sample_id_to_cg_map_;
}

void Merge::ComputeResultEdgeMutations() {
  result_edges_.ReadAll(
      [](auto result_edges, auto& /*result_nodes*/, auto& sample_id_to_cg_map,
         auto result_dag) {
        SeqForEach(result_edges, [&](auto& i) {
          auto& [label, edge_id] = i;
          Assert(label.GetParent().GetCompactGenome());
          const CompactGenome& parent = *label.GetParent().GetCompactGenome();

          // Check NodeLabel's SampleId, not DAG structure, because a node might be
          // structurally a leaf in the result DAG but was an internal node in the source
          if (not label.GetChild().GetSampleId().empty()) {
            const CompactGenome& child =
                sample_id_to_cg_map.at(label.GetChild().GetSampleId().ToString());

            result_dag.Get(edge_id).SetEdgeMutations(CompactGenome::ToEdgeMutations(
                result_dag.GetReferenceSequence(), parent, child));
          } else {
            const CompactGenome& child = *label.GetChild().GetCompactGenome();
            result_dag.Get(edge_id).SetEdgeMutations(CompactGenome::ToEdgeMutations(
                result_dag.GetReferenceSequence(), parent, child));
          }
        });
      },
      GetResultNodes(), SampleIdToCGMap(), ResultDAG());
}

bool Merge::ContainsLeafset(const LeafSet& leafset) const {
  return all_leaf_sets_.find(leafset) != nullptr;
}

template <typename DAGSRange, typename NodeLabelsContainer>
void Merge::MergeCompactGenomes(size_t i, const DAGSRange& dags, NodeId below,
                                std::vector<NodeLabelsContainer>& dags_labels) {
  auto dag = GetFullDAG(dags.at(i));
  dag.AssertUA();
  auto& labels = dags_labels.at(i);
  for (auto node : dag.Const().GetNodes()) {
    labels.insert({node, {}});  // TODO Initialize
  }

  for (auto node : dag.Const().GetNodes()) {
    if (below.value != NoId and node.IsUA()) {
      continue;
    }
    if (node.IsLeaf()) {
      Assert(node.Const().HaveSampleId());
      auto id_iter =
          ResultDAG().template AsFeature<Deduplicate<SampleId>>().AddDeduplicated(
              SampleId::Make(node.Const().GetSampleId().value()));
      labels.at(node).SetSampleId(*id_iter.first);
    } else {
      if (not node.GetCompactGenome().empty()) {
        auto cg_iter = ResultDAG()
                           .template AsFeature<Deduplicate<CompactGenome>>()
                           .AddDeduplicated(node.GetCompactGenome());
        labels.at(node).SetCompactGenome(cg_iter.first);
      }
    }
  }
  for (auto leaf_node : dag.Const().GetLeafs()) {
    Assert(leaf_node.HaveSampleId());
    std::string sid{leaf_node.GetSampleId().value()};
    sample_id_to_cg_map_.insert({sid, leaf_node.GetCompactGenome().Copy(&leaf_node)});
  }
}

template <typename DAGSRange, typename NodeLabelsContainer>
void Merge::ComputeLeafSets(size_t i, const DAGSRange& dags, NodeId below,
                            std::vector<NodeLabelsContainer>& dags_labels) {
  auto dag = GetFullDAG(dags.at(i));
  NodeLabelsContainer& labels = dags_labels.at(i);
  labels.reserve(dag.GetNodesCount());
  using ComputedLSType = IdContainer<NodeId, LeafSet, IdContinuity::Sparse,
                                     Ordering::Ordered>;  // FIXME Dense
  ComputedLSType computed_ls = LeafSet::ComputeLeafSets<ComputedLSType>(dag, labels);
  for (auto node : dag.GetNodes()) {
    if (below.value != NoId and node.IsUA()) {
      continue;
    }
    auto& label = labels.at(node);
    auto* ls = computed_ls.At(node);
    if (ls != nullptr and not ls->empty()) {
      auto ls_iter = all_leaf_sets_.insert(std::move(*ls));
      label.SetLeafSet(std::addressof(ls_iter.first));
      Assert(not label.empty());
    }
  }
}

template <typename DAGSRange, typename NodeLabelsContainer>
void Merge::MergeNodes(size_t i, const DAGSRange& dags, NodeId below,
                       const std::vector<NodeLabelsContainer>& dags_labels,
                       std::atomic<size_t>& node_id) {
  auto&& dag = dags.at(i);
  auto& labels = dags_labels.at({i});
  for (auto node : dag.GetNodes()) {
    if (below.value != NoId and node.IsUA()) {
      continue;
    }
    auto& label = labels.at(node);
    Assert(not label.empty());

    NodeId orig_id = [this, &node_id, &label]() {
      NodeId new_id;
      result_nodes_.insert({label, new_id},
                           [&new_id, &node_id, &label, this](auto&& ins_pair) {
                             if (ins_pair.second) {
                               new_id.value = node_id.fetch_add(1);
                               ins_pair.first = new_id;
                               result_node_labels_.insert_or_assign(new_id, label);
                             } else {
                               new_id.value = ins_pair.first.value;
                             }
                           });
      return new_id;
    }();

    if constexpr (std::remove_reference_t<decltype(dag)>::
                      template contains_element_feature<Component::Node, MappedNodes>) {
      dag.Get(node).SetOriginalId(orig_id);
    }
  }
}

template <typename DAGSRange, typename NodeLabelsContainer>
void Merge::MergeEdges(size_t i, const DAGSRange& dags, NodeId below,
                       const std::vector<NodeLabelsContainer>& dags_labels,
                       Reduction<std::vector<Merge::AddedEdge>>& added_edges) {
  auto&& dag = dags.at(i);
  const NodeLabelsContainer& labels = dags_labels.at(i);
  for (auto edge : dag.GetEdges()) {
    if (below.value != NoId and edge.IsUA()) {
      continue;
    }
    const auto& parent_label = labels.at(edge.GetParentId());
    const auto& child_label = labels.at(edge.GetChildId());
    Assert(not parent_label.empty());
    Assert(not child_label.empty());
    Assert(result_nodes_.find(parent_label) != nullptr);
    Assert(result_nodes_.find(child_label) != nullptr);
    bool ins = result_edges_.insert({{parent_label, child_label}, {}}).second;
    if (ins) {
      added_edges.AddElement(
          [](auto& ae, auto& par_lbl, auto& ch_lbl) {
            ae.push_back({{par_lbl, ch_lbl}, {}, {}, {}, {}});
          },
          parent_label, child_label);
    }
  }
}

void Merge::BuildResult(size_t i, std::vector<Merge::AddedEdge>& added_edges,
                        std::atomic<size_t>& edge_id) {
  auto& [edge, id, parent_id, child_id, clade] = added_edges.at(i);
  id = {edge_id.fetch_add(1)};
  NodeLabel result_parent_label = edge.GetParent();
  NodeLabel result_child_label = edge.GetChild();
  NodeId result_parent_id = result_nodes_.at(result_parent_label);
  NodeId result_child_id = result_nodes_.at(result_child_label);
  Assert(result_parent_id.value != NoId);
  Assert(result_child_id.value != NoId);
  Assert(result_parent_id != result_child_id);
  Assert(result_parent_id < ResultDAG().template GetNextAvailableNodeId<MergeDAG>());
  Assert(result_child_id < ResultDAG().template GetNextAvailableNodeId<MergeDAG>());
  parent_id = result_parent_id;
  child_id = result_child_id;
  clade = edge.ComputeCladeIdx();
  Assert(clade.value != NoId);
  ResultDAG().Get(id).Set(result_parent_id, result_child_id, clade);
  ResultDAG().Get(result_parent_id) = result_parent_label.GetCompactGenome();
  ResultDAG().Get(result_child_id) = result_child_label.GetCompactGenome();
  ResultDAG().Get(result_parent_id) = result_parent_label.GetSampleId();
  ResultDAG().Get(result_child_id) = result_child_label.GetSampleId();
  result_edges_.at(edge) = id;
  ComputeResultEdgeMutations(ResultDAG().Get(id), edge);
}

template <typename Edge>
void Merge::ComputeResultEdgeMutations(Edge edge, const EdgeLabel& label) {
  Assert(label.GetParent().GetCompactGenome());
  const CompactGenome& parent = *label.GetParent().GetCompactGenome();

  if (edge.GetChild().GetCladesCount() < 1) {
    std::string sid = label.GetChild().GetSampleId().ToString();
    const CompactGenome* child = sample_id_to_cg_map_.find(sid);
    if (child == nullptr) {
      return;
    }
    edge.SetEdgeMutations(CompactGenome::ToEdgeMutations(
        ResultDAG().GetReferenceSequence(), parent, *child));
  } else {
    const CompactGenome& child = *label.GetChild().GetCompactGenome();
    edge.SetEdgeMutations(CompactGenome::ToEdgeMutations(
        ResultDAG().GetReferenceSequence(), parent, child));
  }
}

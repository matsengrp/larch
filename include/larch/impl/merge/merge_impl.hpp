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

  const bool was_empty = ResultDAG().IsEmpty();

  std::vector<size_t> idxs;
  idxs.resize(dags.size());
  std::iota(idxs.begin(), idxs.end(), 0);

  ParallelForEach(idxs, [&](size_t i) {
    dags.at(i).GetRoot().Validate(true, not dags.at(i).IsTree());
  });

  std::vector<NodeLabelsContainer> dags_labels;
  dags_labels.resize(dags.size());

  ParallelForEach(idxs, [&](size_t i) {
    MergeCompactGenomes(i, dags, below, dags_labels, ResultDAG());
  });

  ParallelForEach(idxs, [&](size_t i) {
    ComputeLeafSets(i, dags, below, dags_labels, all_leaf_sets_);
  });

  ParallelForEach(idxs, [&](size_t i) {
    for (auto node : dags.at(i).GetNodes()) {
      Assert(not dags_labels.at(i).at(node).Empty());
    }
  });

  std::atomic<size_t> node_id{ResultDAG().GetNextAvailableNodeId().value};
  tbb::concurrent_vector<std::tuple<EdgeLabel, EdgeId, NodeId, NodeId, CladeIdx>>
      added_edges;
  ParallelForEach(idxs, [&](size_t i) {
    MergeNodes(i, dags, below, dags_labels, result_nodes_, result_node_labels_,
               node_id);
  });

  result_nodes_.Read([](auto& result_nodes) {
    for (auto& i : result_nodes) {
      Assert(i.second.value != NoId);
    }
  });

  ParallelForEach(idxs, [&](size_t i) {
    MergeEdges(i, dags, below, dags_labels, result_nodes_, result_edges_, added_edges);
  });

  ResultDAG().InitializeNodes(
      result_nodes_.Read([](auto& result_nodes) { return result_nodes.size(); }));
  std::atomic<size_t> edge_id{ResultDAG().GetNextAvailableEdgeId().value};
  ResultDAG().InitializeEdges(
      result_edges_.Read([](auto& result_edges) { return result_edges.size(); })

  );
  idxs.resize(added_edges.size());
  std::iota(idxs.begin(), idxs.end(), 0);
  ParallelForEach(idxs, [&](size_t i) {
    BuildResult(i, added_edges, edge_id, result_nodes_, result_edges_, ResultDAG());
  });

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

  Assert(result_nodes_.Read([](auto& result_nodes) { return result_nodes.size(); }) ==
         ResultDAG().GetNodesCount());
  Assert(result_node_labels_.Read([](auto& result_node_labels) {
    return result_node_labels.size();
  }) == ResultDAG().GetNodesCount());
  Assert(result_edges_.Read([](auto& result_edges) { return result_edges.size(); }) ==
         ResultDAG().GetEdgesCount());
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

MergeDAG Merge::GetResult() const { return result_dag_storage_.View().Const(); }

MutableMergeDAG Merge::ResultDAG() { return result_dag_storage_.View(); }

const ConcurrentUnorderedMap<NodeLabel, NodeId>& Merge::GetResultNodes() const {
  return result_nodes_;
}

const ConcurrentUnorderedMap<NodeId, NodeLabel>& Merge::GetResultNodeLabels() const {
  return result_node_labels_;
}

void Merge::ComputeResultEdgeMutations() {
  // TODO parallel
  result_edges_.Read(
      [](auto& result_edges, auto result_dag) {
        for (auto& [label, edge_id] : result_edges) {
          Assert(label.GetParent().GetCompactGenome());
          Assert(label.GetChild().GetCompactGenome());
          const CompactGenome& parent = *label.GetParent().GetCompactGenome();
          const CompactGenome& child = *label.GetChild().GetCompactGenome();
          result_dag.Get(edge_id).SetEdgeMutations(CompactGenome::ToEdgeMutations(
              result_dag.GetReferenceSequence(), parent, child));
        }
      },
      ResultDAG());
}

bool Merge::ContainsLeafset(const LeafSet& leafset) const {
  return all_leaf_sets_.find(leafset) != all_leaf_sets_.end();
}

namespace {

template <typename DAG>
auto GetFullDAG(DAG dag) {
  static_assert(DAG::role == Role::View);
  static_assert(DAG::component == Component::DAG);
  return dag;
}

template <typename DAG, template <typename, typename> typename Base>
auto GetFullDAG(DAGView<FragmentStorage<DAG>, Base> dag) {
  return dag.GetStorage().GetTargetStorage().View();
}

template <typename DAG, template <typename, typename> typename Base>
auto GetFullDAG(DAGView<const FragmentStorage<DAG>, Base> dag) {
  return dag.GetStorage().GetTargetStorage().View();
}

}  // namespace

template <typename DAGSRange>
void Merge::MergeCompactGenomes(size_t i, const DAGSRange& dags, NodeId below,
                                std::vector<NodeLabelsContainer>& dags_labels,
                                MutableMergeDAG result_dag) {
  auto dag = GetFullDAG(dags.at(i));
  dag.AssertUA();
  auto& labels = dags_labels.at(i);
  for (auto node : dag.Const().GetNodes()) {
    labels.insert({node, {}});
  }

  for (auto node : dag.Const().GetNodes()) {
    if (below.value != NoId and node.IsUA()) {
      continue;
    }
    if (not node.IsLeaf()) {
      auto cg_iter =
          result_dag.template AsFeature<Deduplicate<CompactGenome>>().AddDeduplicated(
              node.GetCompactGenome());
      labels.at(node).SetCompactGenome(cg_iter.first);
    }
  }
  for (auto leaf : dag.GetLeafs()) {
    Assert(leaf.Const().HaveSampleId());
    auto id_iter =
        result_dag.template AsFeature<Deduplicate<SampleId>>().AddDeduplicated(
            leaf.Const().GetSampleId());
    labels.at(leaf).SetSampleId(id_iter.first);
  }
}

template <typename DAGSRange>
void Merge::ComputeLeafSets(size_t i, const DAGSRange& dags, NodeId below,
                            std::vector<NodeLabelsContainer>& dags_labels,
                            ConcurrentUnorderedSet<LeafSet>& all_leaf_sets) {
  auto dag = GetFullDAG(dags.at(i));
  NodeLabelsContainer& labels = dags_labels.at(i);
  ContiguousMap<NodeId, LeafSet> computed_ls = LeafSet::ComputeLeafSets(dag, labels);
  for (auto node : dag.GetNodes()) {
    if (below.value != NoId and node.IsUA()) {
      continue;
    }
    auto& ls = computed_ls.at(node);
    auto& label = labels.at(node);
    auto ls_iter = all_leaf_sets.insert(std::move(ls));
    label.SetLeafSet(std::addressof(*ls_iter.first));
    Assert(not label.Empty());
  }
}

template <typename DAGSRange>
void Merge::MergeNodes(size_t i, const DAGSRange& dags, NodeId below,
                       std::vector<NodeLabelsContainer>& dags_labels,
                       ConcurrentUnorderedMap<NodeLabel, NodeId>& result_nodes,
                       ConcurrentUnorderedMap<NodeId, NodeLabel>& result_node_labels,
                       std::atomic<size_t>& node_id) {
  auto&& dag = dags.at(i);
  auto& labels = dags_labels.at({i});
  for (auto node : dag.GetNodes()) {
    if (below.value != NoId and node.IsUA()) {
      continue;
    }
    auto& label = labels.at(node);
    Assert(not label.Empty());

    NodeId orig_id = result_nodes.Write(
        [&node_id, &result_node_labels, &label](auto& result_nodes_w) {
          NodeId new_id;
          auto ins_pair = result_nodes_w.insert({label, new_id});
          if (ins_pair.second) {
            new_id.value = node_id.fetch_add(1);
            ins_pair.first->second = new_id;
            result_node_labels.Write([&new_id, &label](auto& result_node_labels_w) {
              result_node_labels_w[new_id] = label;
            });
          } else {
            new_id.value = ins_pair.first->second.value;
          }
          Assert(new_id.value != NoId);
          return new_id;
        });

    if constexpr (std::remove_reference_t<decltype(dag)>::
                      template contains_element_feature<Component::Node, MappedNodes>) {
      dag.Get(node).SetOriginalId(orig_id);
    }
  }
}

template <typename DAGSRange>
void Merge::MergeEdges(
    size_t i, const DAGSRange& dags, NodeId below,
    std::vector<NodeLabelsContainer>& dags_labels,
    const ConcurrentUnorderedMap<NodeLabel, NodeId>& result_nodes,
    ConcurrentUnorderedMap<EdgeLabel, EdgeId>& result_edges,
    tbb::concurrent_vector<std::tuple<EdgeLabel, EdgeId, NodeId, NodeId, CladeIdx>>&
        added_edges) {
  auto&& dag = dags.at(i);
  const NodeLabelsContainer& labels = dags_labels.at(i);
  for (auto edge : dag.GetEdges()) {
    if (below.value != NoId and edge.IsUA()) {
      continue;
    }
    const auto& parent_label = labels.at(edge.GetParentId());
    const auto& child_label = labels.at(edge.GetChildId());
    Assert(not parent_label.Empty());
    Assert(not child_label.Empty());
    result_nodes.Read([&parent_label, &child_label](auto& result_nodes_r) {
      Assert(result_nodes_r.find(parent_label) != result_nodes_r.end());
      Assert(result_nodes_r.find(child_label) != result_nodes_r.end());
    });
    bool ins = result_edges.Write([&parent_label, &child_label](auto& result_edges_w) {
      return result_edges_w.insert({{parent_label, child_label}, {}}).second;
    });
    if (ins) {
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
    ConcurrentUnorderedMap<EdgeLabel, EdgeId>& result_edges,
    MutableMergeDAG result_dag) {
  auto& [edge, id, parent_id, child_id, clade] = added_edges.at(i);
  id = {edge_id.fetch_add(1)};
  auto [result_parent_label, result_parent_id, result_child_label, result_child_id] =
      result_nodes.Read([&edge](auto& result_nodes_r) {
        auto parent = result_nodes_r.find(edge.GetParent());
        auto child = result_nodes_r.find(edge.GetChild());
        Assert(parent != result_nodes_r.end());
        Assert(child != result_nodes_r.end());
        return std::make_tuple(std::ref(parent->first), parent->second,
                               std::ref(child->first), child->second);
      });
  Assert(result_parent_id.value != NoId);
  Assert(result_child_id.value != NoId);
  Assert(result_parent_id != result_child_id);
  if (result_dag.Get(result_child_id).IsLeaf()) {
    Assert(not result_child_label.GetSampleId()->IsEmpty());
  }
  Assert(result_parent_id < result_dag.GetNextAvailableNodeId());
  Assert(result_child_id < result_dag.GetNextAvailableNodeId());
  parent_id = result_parent_id;
  child_id = result_child_id;
  clade = edge.ComputeCladeIdx();
  Assert(clade.value != NoId);
  result_dag.Get(id).Set(result_parent_id, result_child_id, clade);
  result_dag.Get(result_parent_id) = result_parent_label.GetCompactGenome();
  result_dag.Get(result_child_id) = result_child_label.GetCompactGenome();
  result_dag.Get(result_parent_id) = result_parent_label.GetSampleId();
  result_dag.Get(result_child_id) = result_child_label.GetSampleId();
  result_edges.Write([&edge, &id](auto& result_edges_w) {
    auto result_edge_it = result_edges_w.find(edge);
    Assert(result_edge_it != result_edges_w.end());
    result_edge_it->second = id;
  });
}

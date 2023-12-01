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

  ParallelForEach(idxs, [&](size_t i) {
    dags.at(i).GetRoot().Validate(true, not dags.at(i).IsTree());
  });

  constexpr IdContinuity id_continuity = std::remove_reference_t<decltype(dags.at(
      0))>::template id_continuity<Component::Node>;
  std::vector<IdContainer<NodeId, NodeLabel, id_continuity, Ordering::Ordered>>
      dags_labels;
  dags_labels.resize(dags.size());

  ParallelForEach(idxs, [&](size_t i) {
    MergeCompactGenomes(i, dags, below, dags_labels, sample_id_to_cg_map_, ResultDAG());
  });

  ParallelForEach(idxs, [&](size_t i) {
    ComputeLeafSets(i, dags, below, dags_labels, all_leaf_sets_);
  });

#ifndef NDEBUG
  ParallelForEach(idxs, [&](size_t i) {
    for (auto node : dags.at(i).GetNodes()) {
      Assert(not dags_labels.at(i).at(node).empty());
    }
  });
#endif

  std::atomic<size_t> node_id{ResultDAG().GetNextAvailableNodeId().value};
  ParallelForEach(idxs, [&](size_t i) {
    MergeNodes(i, dags, below, dags_labels, result_nodes_, result_node_labels_,
               node_id);
  });

#ifndef NDEBUG
  result_nodes_.Read([](auto& result_nodes) {
    for (auto& i : result_nodes) {
      Assert(i.second.value != NoId);
    }
  });
#endif

  Reduction<std::vector<AddedEdge>> added_edges_reduction;
  ParallelForEach(idxs, [&](size_t i) {
    MergeEdges(i, dags, below, dags_labels, result_nodes_, result_edges_,
               added_edges_reduction);
  });
  std::vector<AddedEdge> added_edges;
  added_edges_reduction.Gather(
      [](auto& ae_threads, auto& result) {
        for (auto& ae : ae_threads) {
          result.insert(result.end(), ae.second.begin(), ae.second.end());
        }
      },
      added_edges);
  added_edges_reduction.clear();

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

const ConcurrentUnorderedMap<NodeLabel, NodeId>& Merge::GetResultNodes() const {
  return result_nodes_;
}

const ConcurrentUnorderedMap<NodeId, NodeLabel>& Merge::GetResultNodeLabels() const {
  return result_node_labels_;
}

const ConcurrentUnorderedMap<std::string, CompactGenome>& Merge::SampleIdToCGMap()
    const {
  return sample_id_to_cg_map_;
}

void Merge::ComputeResultEdgeMutations() {
  // TODO parallel
  result_edges_.Read(
      [](auto& result_edges, auto& result_nodes, auto& sample_id_to_cg_map,
         auto result_dag) {
        for (auto& [label, edge_id] : result_edges) {
          Assert(label.GetParent().GetCompactGenome());
          const CompactGenome& parent = *label.GetParent().GetCompactGenome();

          auto child_node = result_nodes.Read(
              [](auto& result_nodes_r, auto child_label) {
                return result_nodes_r.at(child_label);
              },
              label.GetChild());

          if (result_dag.Get(child_node).IsLeaf()) {
            const CompactGenome& child = sample_id_to_cg_map.Read(
                [](auto& sample_id_to_cg_map_r, auto node_label) {
                  return sample_id_to_cg_map_r.at(node_label.GetSampleId()->ToString())
                      .Copy();
                },
                label.GetChild());

            result_dag.Get(edge_id).SetEdgeMutations(CompactGenome::ToEdgeMutations(
                result_dag.GetReferenceSequence(), parent, child));
          } else {
            const CompactGenome& child = *label.GetChild().GetCompactGenome();
            result_dag.Get(edge_id).SetEdgeMutations(CompactGenome::ToEdgeMutations(
                result_dag.GetReferenceSequence(), parent, child));
          }
        }
      },
      GetResultNodes(), SampleIdToCGMap(), ResultDAG());
}

bool Merge::ContainsLeafset(const LeafSet& leafset) const {
  return all_leaf_sets_.Read(
      [](auto& read, const LeafSet& ls) { return read.find(ls) != read.end(); },
      leafset);
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

template <typename DAGSRange, typename NodeLabelsContainer>
void Merge::MergeCompactGenomes(
    size_t i, const DAGSRange& dags, NodeId below,
    std::vector<NodeLabelsContainer>& dags_labels,
    ConcurrentUnorderedMap<std::string, CompactGenome>& sample_id_to_cg_map,
    MutableMergeDAG result_dag) {
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
          result_dag.template AsFeature<Deduplicate<SampleId>>().AddDeduplicated(
              node.Const().GetSampleId());
      labels.at(node).SetSampleId(id_iter.first);
    } else {
      if (not node.GetCompactGenome().empty()) {
        auto cg_iter =
            result_dag.template AsFeature<Deduplicate<CompactGenome>>().AddDeduplicated(
                node.GetCompactGenome());
        labels.at(node).SetCompactGenome(cg_iter.first);
      }
    }
  }
  for (auto leaf_node : dag.Const().GetLeafs()) {
    Assert(leaf_node.HaveSampleId());
    std::string sid = leaf_node.GetSampleId().value();
    sample_id_to_cg_map.Write([&sid, &leaf_node](auto& sample_id_to_cg_map_w) {
      sample_id_to_cg_map_w.insert({sid, leaf_node.GetCompactGenome().Copy()});
    });
  }
}

template <typename DAGSRange, typename NodeLabelsContainer>
void Merge::ComputeLeafSets(size_t i, const DAGSRange& dags, NodeId below,
                            std::vector<NodeLabelsContainer>& dags_labels,
                            ConcurrentUnorderedSet<LeafSet>& all_leaf_sets) {
  auto dag = GetFullDAG(dags.at(i));
  NodeLabelsContainer& labels = dags_labels.at(i);
  using ComputedLSType =
      IdContainer<NodeId, LeafSet, IdContinuity::Sparse, Ordering::Ordered>;
  ComputedLSType computed_ls = LeafSet::ComputeLeafSets<ComputedLSType>(dag, labels);
  for (auto node : dag.GetNodes()) {
    if (below.value != NoId and node.IsUA()) {
      continue;
    }
    auto& label = labels.at(node);
    auto& ls = computed_ls.at(node);
    if (not ls.empty()) {
      all_leaf_sets.Write(
          [](auto& write, auto& lset, auto& lbl) {
            auto ls_iter = write.insert(std::move(lset));
            lbl.SetLeafSet(std::addressof(*ls_iter.first));
          },
          ls, label);
    }
    Assert(not label.empty());
  }
}

template <typename DAGSRange, typename NodeLabelsContainer>
void Merge::MergeNodes(size_t i, const DAGSRange& dags, NodeId below,
                       const std::vector<NodeLabelsContainer>& dags_labels,
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
    Assert(not label.empty());

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

template <typename DAGSRange, typename NodeLabelsContainer>
void Merge::MergeEdges(
    size_t i, const DAGSRange& dags, NodeId below,
    const std::vector<NodeLabelsContainer>& dags_labels,
    [[maybe_unused]] const ConcurrentUnorderedMap<NodeLabel, NodeId>& result_nodes,
    ConcurrentUnorderedMap<EdgeLabel, EdgeId>& result_edges,
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
#ifndef NDEBUG
    result_nodes.Read([&parent_label, &child_label](auto& result_nodes_r) {
      Assert(result_nodes_r.find(parent_label) != result_nodes_r.end());
      Assert(result_nodes_r.find(child_label) != result_nodes_r.end());
    });
#endif
    bool ins = result_edges.Write([&parent_label, &child_label](auto& result_edges_w) {
      return result_edges_w.insert({{parent_label, child_label}, {}}).second;
    });
    if (ins) {
      added_edges.Add(
          [](auto& ae, auto& par_lbl, auto& ch_lbl) {
            ae.push_back({{par_lbl, ch_lbl}, {}, {}, {}, {}});
          },
          parent_label, child_label);
    }
  }
}

void Merge::BuildResult(size_t i, std::vector<Merge::AddedEdge>& added_edges,
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

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsMATRoot() const {
  auto& node = static_cast<const CRTP&>(*this);
  auto old = node.GetOld();
  if (old.IsUA()) {
    return false;
  }
  return old.GetSingleParent().GetParent().IsUA();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsMoveSource() const {
  auto& node = static_cast<const CRTP&>(*this);
  auto move_src = node.GetDAG().GetMoveSource();
  return (move_src.GetId() == node.GetId());
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsMoveTarget() const {
  auto& node = static_cast<const CRTP&>(*this);
  auto move_dst = node.GetDAG().GetMoveTarget();
  return (move_dst.GetId() == node.GetId());
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsMoveNew() const {
  auto& node = static_cast<const CRTP&>(*this);
  return node.GetDAG().GetMoveNew().GetId() == node.GetId();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsLCAAncestor() const {
  auto& node = static_cast<const CRTP&>(*this);
  auto& ancestors = node.GetDAG().GetLCAAncestors();
  return ancestors.find(node) != ancestors.end();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::HasChangedTopology() const {
  return GetFeatureStorage(this).get().has_changed_topology_;
}

template <typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalNode, CRTP, Tag>::GetOld() const {
  auto& node = static_cast<const CRTP&>(*this);
  return node.GetDAG().GetOriginal().Get(node.GetId());
}

template <typename CRTP, typename Tag>
const ContiguousSet<MutationPosition>&
FeatureConstView<HypotheticalNode, CRTP, Tag>::GetChangedBaseSites() const {
  return GetFeatureStorage(this).get().changed_base_sites_;
}

template <typename CRTP, typename Tag>
ContiguousSet<MutationPosition>
FeatureConstView<HypotheticalNode, CRTP, Tag>::GetSitesWithChangedFitchSets() const {
  auto& node = static_cast<const CRTP&>(*this);
  ContiguousSet<MutationPosition> result;
  std::optional<ContiguousMap<MutationPosition, Mutation_Count_Change>> fitch_set_map =
      node.GetFitchSetParts().second;
  if (fitch_set_map.has_value()) {
    for (auto map_pair : fitch_set_map.value()) {
      result.insert(map_pair.first);
    }
  }
  return result;
}

template <typename CRTP, typename Tag>
std::pair<MAT::Mutations_Collection,
          std::optional<ContiguousMap<MutationPosition, Mutation_Count_Change>>>
FeatureConstView<HypotheticalNode, CRTP, Tag>::GetFitchSetParts() const {
  auto& node = static_cast<const CRTP&>(*this);
  auto dag = node.GetDAG();
  return dag.GetBackend().GetFitchSetParts(dag, node.GetId(), node.GetOld().IsLeaf(),
                                           node.IsMoveTarget(), node.IsMoveNew());
}

template <typename CRTP, typename Tag>
FitchSet FeatureConstView<HypotheticalNode, CRTP, Tag>::GetFitchSetAtSite(
    MutationPosition site) const {
  auto node = static_cast<const CRTP&>(*this).Const();
  auto dag = node.GetDAG();
  return dag.GetBackend().GetFitchSetAtSite(dag, node.GetId(), site,
                                            node.GetOld().IsLeaf(), node.IsMoveTarget(),
                                            node.IsMoveNew());
}

template <typename CRTP, typename Tag>
ContiguousSet<MutationPosition>
FeatureConstView<HypotheticalNode, CRTP, Tag>::GetParentChangedBaseSites() const {
  auto& node = static_cast<const CRTP&>(*this);
  auto dag = node.GetDAG();
  if (node.IsMoveSource()) {
    // this node used to be below the old source parent, so we need
    // differences between the new parent's new cg, and the old cg of the old
    // parent of the source node.
    const CompactGenome& old_parent_cg = dag.GetOldSourceParent().GetCompactGenome();
    // Parent must be the new node:
    const CompactGenome& new_parent_cg =
        node.GetSingleParent().GetParent().GetCompactGenome();
    // Imaginary method DifferingSites returns sites at which new_parent_cg
    // and old_parent_cg don't have the same base.
    return old_parent_cg.DifferingSites(new_parent_cg);
  } else if (node.IsMoveTarget()) {
    // if this node is the target node, then the old parent of this node is now its
    // grandparent, so we want to check which sites differ between the old parent and
    // the new parent (which is the new node).
    const CompactGenome& old_parent_cg =
        dag.GetMoveTarget().GetOld().GetSingleParent().GetParent().GetCompactGenome();
    const CompactGenome& new_parent_cg =
        node.GetSingleParent().GetParent().GetCompactGenome();
    // Imaginary method DifferingSites returns sites at which new_parent_cg
    // and old_parent_cg don't have the same base.
    return old_parent_cg.DifferingSites(new_parent_cg);
  } else {
    return node.GetSingleParent().GetParent().GetChangedBaseSites().Copy();
  }
}

template <typename CRTP, typename Tag>
CompactGenome FeatureConstView<HypotheticalNode, CRTP, Tag>::ComputeNewCompactGenome()
    const {
  auto node = static_cast<const CRTP&>(*this).Const();
  ContiguousSet<MutationPosition> changed_base_sites =
      node.GetSitesWithChangedFitchSets();
  const CompactGenome& old_cg = node.GetOld().GetCompactGenome();
  ContiguousMap<MutationPosition, MutationBase> cg_changes;
  if (node.IsMATRoot()) {
    // If this node is the root node of the tree, we don't have to worry
    // about parent bases, we just
    // choose a base from each changed fitch set, preferring the base that was
    // already in that site before the SPR move
    ContiguousSet<MutationPosition> focus_sites = GetSitesWithChangedFitchSets();
    for (auto site : focus_sites) {
      FitchSet site_fitch_set = GetFitchSetAtSite(site);
      auto oldbase = old_cg.GetBase(site, node.GetDAG().GetReferenceSequence());
      if (not site_fitch_set.find(oldbase.ToChar())) {
        cg_changes.insert({site, site_fitch_set.at(0)});
        changed_base_sites.insert(site);
      }
    }
  } else {
    // If this node is not the root node, we do need to prefer the base
    // from the new compact genome of the parent node. If it's not in the
    // fitch set, then we're free to choose any base in the fitch set, so we
    // choose the base in this node's old compact genome when that base is in
    // the fitch set.
    ContiguousSet<MutationPosition> focus_sites = node.GetSitesWithChangedFitchSets();
    focus_sites.Union(node.GetParentChangedBaseSites());
    for (auto site : focus_sites) {
      FitchSet site_fitch_set = node.GetFitchSetAtSite(site);
      auto oldbase = old_cg.GetBase(site, node.GetDAG().GetReferenceSequence());
      auto parent_base = node.GetSingleParent().GetParent().GetCompactGenome().GetBase(
          site, node.GetDAG().GetReferenceSequence());
      if (site_fitch_set.find(parent_base.ToChar())) {
        if (oldbase != parent_base) {
          cg_changes.insert({site, parent_base});
          changed_base_sites.insert(site);
        }
      } else if (not site_fitch_set.find(oldbase.ToChar())) {
        cg_changes.insert({site, site_fitch_set.at(0)});
        changed_base_sites.insert(site);
      }
    }
  }
  CompactGenome result = old_cg.Copy(static_cast<const CRTP*>(this));
  result.ApplyChanges(cg_changes);
  return result;
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsNonrootAnchorNode() const {
  auto node = static_cast<const CRTP&>(*this).Const();
  if (node.IsMoveSource()) {
    return false;
  }
  if (node.IsMoveTarget()) {
    return false;
  }
  if (node.IsMoveNew() or node.HasChangedTopology()) {
    return false;
  }
  if (not IsLCAAncestor()) {
    return node.GetOld().GetCompactGenome() == node.GetCompactGenome();
  } else {
    return false;
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<HypotheticalNode, CRTP, Tag>::SetHasChangedTopology() const {
  GetFeatureStorage(this).get().has_changed_topology_ = true;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<HypotheticalNode, CRTP, Tag>::PreorderComputeCompactGenome(
    std::vector<NodeId>& result_nodes, std::vector<EdgeId>& result_edges) const {
  auto& node = static_cast<const CRTP&>(*this);
  if (not node.IsUA() and not node.IsMoveNew() and not node.IsLeaf()) {
    if (not node.template IsOverlaid<Deduplicate<CompactGenome>>()) {
      node.template SetOverlay<Deduplicate<CompactGenome>>();
    }
    node = node.ComputeNewCompactGenome();
  }
  result_nodes.push_back(node);
  // If we've reached an anchor node, there's no need to continue down this
  // branch.
  if (node.IsUA() or (not(node.IsNonrootAnchorNode() or result_nodes.size() < 2))) {
    for (auto child : node.GetChildren()) {
      result_edges.push_back(child);
      child.GetChild().PreorderComputeCompactGenome(result_nodes, result_edges);
    }
  }
}

// HypotheticalTree FeatureConstView implementations

template <typename DAG, typename Backend, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::GetMoveLCA() const {
  auto& self = GetFeatureStorage(this).get();
  Assert(self.data_);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.Get(self.data_->backend_.GetMoveLCA());
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::GetMoveSource()
    const {
  auto& self = GetFeatureStorage(this).get();
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.Get(self.data_->backend_.GetMoveSource());
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::GetMoveTarget()
    const {
  auto& self = GetFeatureStorage(this).get();
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.Get(self.data_->backend_.GetMoveTarget());
}

// TODO_DR: Used by larch_usher.cpp
/* CONDENSING CODE: can remove this function (replace all calls to it with calls to
 * `GetMoveSource' */
template <typename DAG, typename Backend, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::GetMoveSources()
    const {
  auto& self = GetFeatureStorage(this).get();
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetUncondensedNodeFromMAT(self.data_->move_.src);
}

// TODO_DR: Used by larch_usher.cpp
/* CONDENSING CODE: can remove this function (replace all calls to it with calls to
 * `GetMoveTarget' */
template <typename DAG, typename Backend, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::GetMoveTargets()
    const {
  auto& self = GetFeatureStorage(this).get();
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetUncondensedNodeFromMAT(self.data_->move_.dst);
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::GetMoveNew() const {
  auto& self = GetFeatureStorage(this).get();
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.Get(self.data_->new_node_);
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP,
                      Tag>::HasUnifurcationAfterMove() const {
  auto& self = GetFeatureStorage(this).get();
  return self.data_->has_unifurcation_after_move_;
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::GetOldSourceParent()
    const {
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetOriginal()
      .Get(dag.GetMoveSource().GetId())
      .GetSingleParent()
      .GetParent();
}

template <typename Node>
static inline bool IsChanged(Node node) {
  return node.GetCompactGenome() != node.GetOld().GetCompactGenome();
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::GetOldestChangedNode()
    const {
  auto& dag = static_cast<const CRTP&>(*this);

  NodeId oldest = [dag]() {
    auto lca = dag.GetOriginal().GetMoveLCA();
    if (lca.GetCladesCount() != 2) {
      return lca;
    }
    if (not lca.ContainsChild(dag.GetMoveSource())) {
      return lca;
    }
    for (auto node : lca.GetChildren() | Transform::GetChild()) {
      auto src = dag.GetMoveSource();
      if (src == node.GetId()) {
        return node;
      }
    }
    Fail("oldest changed node Unreachable");
  }();

  while (not dag.Get(oldest).IsUA() and IsChanged(dag.Get(oldest).Const())) {
    oldest = dag.Get(oldest).GetSingleParent().GetParent();
  }
  return dag.Get(oldest);
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::MakeFragment() const {
  auto& dag = static_cast<const CRTP&>(*this);

  std::vector<NodeId> result_nodes;
  std::vector<EdgeId> result_edges;
  auto oldest_changed = dag.GetOldestChangedNode().GetSingleParent().GetParent();
  if (oldest_changed.IsUA()) {
    // we need to add the UA node as the root anchor node of the fragment,
    // somehow
    // TODO how is this being done, then? (Things seem to be working...)
  } else {
    result_nodes.push_back(oldest_changed.GetSingleParent().GetParent());
    result_edges.push_back(oldest_changed.GetSingleParent());
  }
  oldest_changed.PreorderComputeCompactGenome(result_nodes, result_edges);

  auto collapsed = dag.CollapseEmptyFragmentEdges(result_nodes, result_edges);
  NodeId oldest_node = collapsed.first.front();

#ifdef KEEP_ASSERTS
  for (auto node_id : collapsed.first) {
    auto node = dag.Get(node_id);
    if (not node.IsUA() and not node.IsMoveNew()) {
      if (not node.GetOld().HaveSampleId()) {
        Assert(node.GetCladesCount() > 0);
        for (auto clade : node.GetClades()) {
          Assert(not clade.empty());
        }
      }
    }
  }
#endif

  for (auto node_id : collapsed.first) {
    if (dag.Get(node_id).IsUA()) {
      oldest_node = node_id;
    } else {
      auto parent_edge = dag.Get(node_id).GetSingleParent().GetId();
      if (std::find(collapsed.second.begin(), collapsed.second.end(), parent_edge) ==
          collapsed.second.end()) {
        oldest_node = node_id;
        break;
      }
    }
  }
  // make sure that all children of the most ancestral node are included,
  // not just the child that is the fragment's root.
  for (auto child : dag.Get(oldest_node).GetChildren()) {
    collapsed.first.push_back(child.GetChild());
    collapsed.second.push_back(child);
  }
  collapsed.first |= ranges::actions::sort(std::less<NodeId>{}) |
                     ranges::actions::unique(std::equal_to<NodeId>{});
  collapsed.second |= ranges::actions::sort(std::less<EdgeId>{}) |
                      ranges::actions::unique(std::equal_to<EdgeId>{});
  return AddFragmentStorage(dag, std::move(collapsed.first),
                            std::move(collapsed.second), oldest_node);
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP,
                      Tag>::MakeUncollapsedFragment() const {
  auto& dag = static_cast<const CRTP&>(*this);
  std::vector<NodeId> result_nodes;
  std::vector<EdgeId> result_edges;
  auto oldest_changed = dag.GetOldestChangedNode().GetSingleParent().GetParent();
  NodeId oldest_node = dag.GetOldestChangedNode().GetSingleParent().GetParent().GetId();
  if (oldest_changed.IsUA()) {
    // we need to add the UA node as the root anchor node of the fragment,
    // somehow
    // TODO how is this being done, then? (Things seem to be working...)
  } else {
    result_nodes.push_back(oldest_changed.GetSingleParent().GetParent());
    result_edges.push_back(oldest_changed.GetSingleParent());
    oldest_node = oldest_changed.GetSingleParent().GetParent().GetId();
  }
  oldest_changed.PreorderComputeCompactGenome(result_nodes, result_edges);
  // make sure that all children of the most ancestral node are included,
  // not just the child that is the fragment's root.
  for (auto child : dag.Get(oldest_node).GetChildren()) {
    result_nodes.push_back(child.GetChild());
    result_edges.push_back(child);
  }
  result_nodes |= ranges::actions::sort(std::less<NodeId>{}) |
                  ranges::actions::unique(std::equal_to<NodeId>{});
  result_edges |= ranges::actions::sort(std::less<EdgeId>{}) |
                  ranges::actions::unique(std::equal_to<EdgeId>{});
  return AddFragmentStorage(dag, std::move(result_nodes), std::move(result_edges),
                            oldest_node);
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
std::pair<std::vector<NodeId>, std::vector<EdgeId>>
FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::CollapseEmptyFragmentEdges(
    const std::vector<NodeId>& fragment_nodes,
    const std::vector<EdgeId>& fragment_edges) const {
  auto& dag = static_cast<const CRTP&>(*this);

  /*/
  // BEGIN WORKAROUND CODE FOR AVOIDING COLLAPSING
  IdContainer<NodeId, bool, IdContinuity::Sparse, Ordering::Unordered>
      is_parent_of_collapsible_edge;
  IdContainer<NodeId, bool, IdContinuity::Sparse, Ordering::Unordered>
      is_child_of_collapsible_edge;
  IdContainer<NodeId, bool, IdContinuity::Sparse, Ordering::Unordered>
      node_already_added;
  IdContainer<EdgeId, bool, IdContinuity::Sparse, Ordering::Unordered>
      edge_already_added;
  for (auto e : fragment_edges) {
    auto edge = dag.Get(e);
    auto parent = dag.Get(e).GetParent();
    auto child = dag.Get(e).GetChild();
    auto clade = edge.GetClade().value;
    if (not child.template IsOverlaid<HypotheticalNode>()) {
      child.template SetOverlay<HypotheticalNode>();
    }
    if (not child.template IsOverlaid<Neighbors>()) {
      child.template SetOverlay<Neighbors>();
    }
    if (not parent.template IsOverlaid<HypotheticalNode>()) {
      parent.template SetOverlay<HypotheticalNode>();
    }
    if (not parent.template IsOverlaid<Neighbors>()) {
      parent.template SetOverlay<Neighbors>();
    }
    if (not edge.template IsOverlaid<Endpoints>()) {
      edge.template SetOverlay<Endpoints>();
      if (parent.IsUA()) {
        edge.Set(parent, child, {0});
        child.SetSingleParent(edge);
      } else {
        edge.Set(parent, child, {clade});
        child.SetSingleParent(edge);
      }
    }
    node_already_added.insert_or_assign(parent, true);
    node_already_added.insert_or_assign(child, true);
    edge_already_added.insert_or_assign(edge, true);
  }
  // END WORKAROUND CODE FOR AVOIDING COLLAPSING

  /*/
  // keep track of edges/nodes that are collapsible
  IdContainer<NodeId, bool, IdContinuity::Sparse, Ordering::Unordered>
      is_parent_of_collapsible_edge;
  IdContainer<NodeId, bool, IdContinuity::Sparse, Ordering::Unordered>
      is_child_of_collapsible_edge;
  IdContainer<NodeId, bool, IdContinuity::Sparse, Ordering::Unordered>
      parent_is_in_fragment;
  IdContainer<NodeId, size_t, IdContinuity::Sparse, Ordering::Unordered> clades_count;
  IdContainer<NodeId, bool, IdContinuity::Sparse, Ordering::Unordered>
      node_already_added;
  IdContainer<EdgeId, bool, IdContinuity::Sparse, Ordering::Unordered>
      edge_already_added;
  IdContainer<NodeId, bool, IdContinuity::Sparse, Ordering::Unordered> visited;
  for (auto edge_id : fragment_edges) {
    auto edge = dag.Get(edge_id);
    auto parent = edge.GetParent();
    auto child = edge.GetChild();
    parent_is_in_fragment.insert({child, true});
    if (not(parent.Const().GetCompactGenome() != child.Const().GetCompactGenome() or
            child.IsLeaf() or parent.IsUA() or child.IsNonrootAnchorNode() or
            parent.IsMoveNew() or child.IsMoveNew())) {
      is_parent_of_collapsible_edge.insert({parent, true});
      is_child_of_collapsible_edge.insert({child, true});
    }
  }
  auto FindFinalParent = [&](auto&& self, NodeId node_id) -> std::pair<EdgeId, NodeId> {
    if (dag.Get(node_id) == dag.GetRoot()) {
      return {EdgeId{NoId}, node_id};
    }
    auto parent_edge = dag.Get(node_id).GetSingleParent();
    auto parent_node = parent_edge.GetParent();
    if (is_child_of_collapsible_edge[node_id]) {
      return self(self, parent_node);
    } else if (is_child_of_collapsible_edge[parent_node]) {
      auto grandparent = self(self, parent_node);
      return {grandparent.first, dag.Get(grandparent.first).GetChild()};
    } else {
      return {parent_edge, parent_node};
    }
  };
  // visit each of the fragment's nodes, and, if it is not a collapsible node, add it
  // to the set of nodes that are in the collapsed segment. A node is collapsible if
  // it is the child of a collapsible edge (i.e. an edge with no mutations on it)
  for (auto node_id : fragment_nodes) {
    auto this_node = dag.Get(node_id);
    if ((this_node.IsNonrootAnchorNode() or this_node.IsLeaf() or
         (not is_child_of_collapsible_edge[node_id]))) {
      node_already_added.insert_or_assign(node_id, true);
      if (is_parent_of_collapsible_edge[node_id] and
          not parent_is_in_fragment[node_id]) {
        auto grandparent_edge = this_node.GetSingleParent();

        node_already_added.insert_or_assign(grandparent_edge.GetParent(), true);
        edge_already_added.insert_or_assign(grandparent_edge, true);

        for (auto child_edge : grandparent_edge.GetParent().GetChildren()) {
          if (child_edge != grandparent_edge) {
            if (not node_already_added[child_edge.GetChild().GetId()]) {
              node_already_added.insert_or_assign(child_edge.GetChild(), true);
            }
            if (not edge_already_added[child_edge.GetId()]) {
              edge_already_added.insert_or_assign(child_edge, true);
            }
          }
        }
      }
      if (parent_is_in_fragment[node_id]) {
        auto parent_pair = FindFinalParent(FindFinalParent, node_id);
        auto parent_edge = dag.Get(node_id).GetSingleParent();
        auto parent_node = dag.Get(parent_pair.second);
        auto old_parent_node = this_node.GetSingleParent().GetParent();
        // if this node's parent is changed due to collapsing
        if (parent_node != old_parent_node) {
          auto gp_edge_id = parent_node.GetParentsCount() < 1
                                ? EdgeId{NoId}
                                : parent_node.GetSingleParent();

          std::vector<EdgeId> sibling_edges;
          for (auto sib_edge : parent_node.GetChildren()) {
            if ((edge_already_added[sib_edge] or
                 not is_child_of_collapsible_edge[sib_edge.GetChild()]) and
                (sib_edge.GetChild() != this_node)) {
              sibling_edges.push_back(sib_edge);
            }
          }

          if (not this_node.template IsOverlaid<HypotheticalNode>()) {
            this_node.template SetOverlay<HypotheticalNode>();
          }
          if (not parent_edge.template IsOverlaid<Endpoints>()) {
            parent_edge.template SetOverlay<Endpoints>();
          }
          if (not this_node.template IsOverlaid<Neighbors>()) {
            this_node.template SetOverlay<Neighbors>();
          }
          if (not parent_node.template IsOverlaid<HypotheticalNode>()) {
            parent_node.template SetOverlay<HypotheticalNode>();
          }
          if (not parent_node.template IsOverlaid<Neighbors>()) {
            parent_node.template SetOverlay<Neighbors>();
          }
          auto clade_ins = clades_count.insert({parent_node, 0});
          size_t clade = clade_ins.first->second;
          bool cleared_parent = false;
          if (clade_ins.second) {
            parent_node.ClearConnections();
            cleared_parent = true;
            if (gp_edge_id.value != NoId) {
              auto gp_edge = dag.Get(gp_edge_id);
              auto gp_edge_parent = gp_edge.GetParent();
              auto gp_edge_clade = gp_edge.GetClade().value;
              if (not gp_edge.template IsOverlaid<Endpoints>()) {
                gp_edge.template SetOverlay<Endpoints>();
              }
              gp_edge.Set(gp_edge_parent, parent_node, {gp_edge_clade});
              parent_node.SetSingleParent(gp_edge);
              if (not node_already_added[gp_edge_parent]) {
                node_already_added.insert_or_assign(gp_edge_parent, true);
              }
              if (not edge_already_added[gp_edge_id]) {
                edge_already_added.insert_or_assign(gp_edge_id, true);
              }
            }
            clade = 0;
            for (auto sib_edge : sibling_edges) {
              auto sib_node = dag.Get(sib_edge).GetChild();
              if (not dag.Get(sib_edge).template IsOverlaid<Endpoints>()) {
                dag.Get(sib_edge).template SetOverlay<Endpoints>();
              }
              if (not(node_already_added[sib_node.GetId()] or
                      edge_already_added[sib_edge])) {
                node_already_added.insert({sib_node, true});
              }
              dag.Get(sib_edge).Set(parent_node, sib_node, {clade});
              parent_node.AddEdge({clade++}, sib_edge, true);
              edge_already_added.insert_or_assign(sib_edge, true);
            }
          }
          Assert(parent_edge.template IsOverlaid<Endpoints>());
          parent_edge.Set(parent_node, this_node, {clade});
          Assert(parent_edge.GetParent().GetId() == parent_node.GetId());
          if ((not edge_already_added[parent_edge.GetId()]) or cleared_parent) {
            parent_node.AddEdge({clade}, parent_edge, true);
            clades_count.insert_or_assign(parent_node, clade + 1);
          }
          this_node.SetSingleParent(parent_edge);
          edge_already_added.insert_or_assign(parent_edge, true);
        } else {
          if (not visited[this_node]) {
            node_already_added.insert_or_assign(this_node, true);
            node_already_added.insert_or_assign(parent_node, true);
            edge_already_added.insert_or_assign(parent_edge, true);
            if (not parent_node.template IsOverlaid<HypotheticalNode>()) {
              parent_node.template SetOverlay<HypotheticalNode>();
            }
            if (not this_node.template IsOverlaid<HypotheticalNode>()) {
              this_node.template SetOverlay<HypotheticalNode>();
            }
            if (not parent_node.template IsOverlaid<Neighbors>()) {
              parent_node.template SetOverlay<Neighbors>();
            }
            if (not this_node.template IsOverlaid<Neighbors>()) {
              this_node.template SetOverlay<Neighbors>();
            }
            visited.insert_or_assign(this_node, true);
          }
        }
      }
    }
  }
  //*/
  std::vector<NodeId> current_nodes;
  std::vector<EdgeId> current_edges;
  for (const auto& nodeadded : node_already_added) {
    current_nodes.push_back(nodeadded.first);
    auto this_node = dag.Get(nodeadded.first);
    for (auto c : this_node.GetChildren()) {
      if (not c.template IsOverlaid<Endpoints>()) {
        auto c_child = c.GetChild();
        auto c_clade = c.GetClade().value;
        c.template SetOverlay<Endpoints>();
        c.Set(this_node, c_child, {c_clade});
      }
      edge_already_added.insert_or_assign(c, true);
    }
  }
  for (const auto& edgeadded : edge_already_added) {
    current_edges.push_back(edgeadded.first);
  }
#ifdef KEEP_ASSERTS
  for (auto node_id : current_nodes) {
    auto node = dag.Get(node_id);
    Assert(node_id.value != NoId);
    if (node != dag.GetRoot()) {
      Assert(node.GetParentsCount() == 1);
    }
  }
  for (auto edge_id : current_edges) {
    if (std::find(current_nodes.begin(), current_nodes.end(),
                  dag.Get(edge_id).GetChild()) != current_nodes.end()) {
      auto parent = dag.Get(edge_id).GetParent();
      auto child = dag.Get(edge_id).GetChild();
      if (not child.template IsOverlaid<HypotheticalNode>()) {
        child.template SetOverlay<HypotheticalNode>();
      }
      if (not parent.IsUA()) {
        Assert(parent.GetParentsCount() == 1);
      }
      if (not parent.template IsOverlaid<HypotheticalNode>()) {
        parent.template SetOverlay<HypotheticalNode>();
      }
      Assert(child.GetParentsCount() == 1);
      Assert(edge_id.value != NoId);
      Assert(std::find(current_nodes.begin(), current_nodes.end(), parent) !=
             current_nodes.end());
      Assert(parent.ContainsChild(child));
      Assert(parent.GetCladesCount() > 0);
    }
  }
#endif
  // Assert(current_nodes.size() == current_edges.size() + 1);

  for (auto node_id : fragment_nodes) {
    if (is_parent_of_collapsible_edge[node_id] and
        is_child_of_collapsible_edge[node_id]) {
      // TODO delete / clear this node
    }
  }

  for (auto edge_id : fragment_edges) {
    if (not edge_already_added[edge_id]) {
      // TODO delete / clear this edge
    }
  }
  return {current_nodes, current_edges};
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
const Backend& FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::GetBackend()
    const {
  auto& self = GetFeatureStorage(this).get();
  Assert(self.data_);
  return self.data_->backend_;
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
Score FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::GetScoreChange()
    const {
  auto& self = GetFeatureStorage(this).get();
  Assert(self.data_);
  return self.data_->backend_.GetScoreChange();
}

// Static empty map for non-matOptimize backends
namespace detail {
inline const ContiguousMap<MATNodePtr,
                           ContiguousMap<MutationPosition, Mutation_Count_Change>>&
GetEmptyFitchSetMap() {
  static ContiguousMap<MATNodePtr,
                       ContiguousMap<MutationPosition, Mutation_Count_Change>>
      empty_map;
  return empty_map;
}
}  // namespace detail

template <typename DAG, typename Backend, typename CRTP, typename Tag>
const ContiguousMap<MATNodePtr, ContiguousMap<MutationPosition, Mutation_Count_Change>>&
FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::GetChangedFitchSetMap()
    const {
  auto& self = GetFeatureStorage(this).get();
  Assert(self.data_);
  if constexpr (std::is_same_v<Backend, MatOptimizeScoringBackend<DAG>>) {
    // For matOptimize backend, return the internal map
    // Note: we need to convert from NodeId-keyed to MATNodePtr-keyed
    // For now, we keep both in the backend
    return detail::GetEmptyFitchSetMap();  // TODO: expose mat_keyed_fitch_set_map_ from
                                           // backend
  } else {
    return detail::GetEmptyFitchSetMap();
  }
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
const ContiguousSet<NodeId>&
FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::GetLCAAncestors() const {
  auto& self = GetFeatureStorage(this).get();
  Assert(self.data_);
  return self.data_->lca_ancestors_;
}

namespace {

template <typename DAGType>
std::pair<NodeId, bool> ApplyMoveImpl(DAGType dag, NodeId lca, NodeId& src,
                                      NodeId& dst) {
  for (auto node : dag.GetNodes()) {
    if (not node.IsUA()) {
      if (not node.IsMATRoot()) {
        Assert(node.GetParentsCount() == 1);
      }
    }
  }
  Assert(not dag.HaveOverlays());
  Assert(dag.IsTree());

  auto first_src_node = dag.Get(src);
  auto first_dst_node = dag.Get(dst);
  auto src_parent_node = first_src_node.GetSingleParent().GetParent();
  auto dst_parent_node = first_dst_node.GetSingleParent().GetParent();
  auto src_parent_edge = first_src_node.GetSingleParent();
  auto dst_parent_edge = first_dst_node.GetSingleParent();
  if (first_src_node.IsTreeRoot() or first_src_node.GetId() == first_dst_node.GetId() or
      first_dst_node.GetSingleParent().GetParent().IsUA()) {
    // no-op
    return {};
  }

  const bool is_sibling_move = src_parent_node.GetId() == dst_parent_node.GetId();
  size_t src_size = 1;
  size_t dst_size = 1;
  const bool has_unifurcation_after_move =
      is_sibling_move ? src_parent_node.GetCladesCount() <= (src_size + dst_size)
                      : src_parent_node.GetCladesCount() <= (src_size + 1);
  if ((is_sibling_move and
       (src_parent_node.GetCladesCount() == (src_size + dst_size))) or
      src_parent_node.IsTreeRoot() or src_parent_node.IsUA() or
      (not is_sibling_move and src_parent_node == lca)) {
    // no-op
    return {};
  }
  if (has_unifurcation_after_move and (first_dst_node == src_parent_node)) {
    return {};
  }
  Assert(src_parent_node.GetCladesCount() > 0);
  Assert(dst_parent_node.GetCladesCount() > 0);

#ifdef KEEP_ASSERTS
  const size_t old_num_nodes = dag.GetNodesCount();
  const size_t old_num_edges = dag.GetEdgesCount();
#endif

  auto new_node = has_unifurcation_after_move ? src_parent_node : dag.AppendNode();
  auto new_edge = has_unifurcation_after_move ? src_parent_node.GetSingleParent()
                                              : dag.AppendEdge();
  auto src_edge = EdgeId{NoId};
  std::vector<EdgeId> src_sibling_edges;
  for (auto e : src_parent_node.GetChildren()) {
    if (src != e.GetChild().GetId()) {
      src_sibling_edges.push_back(e.GetId());
    } else {
      src_edge = e.GetId();
    }
  }
  auto dst_edge = EdgeId{NoId};
  std::vector<EdgeId> dst_sibling_edges;
  for (auto e : dst_parent_node.GetChildren()) {
    if (dst != e.GetChild().GetId()) {
      dst_sibling_edges.push_back(e.GetId());
    } else {
      dst_edge = e.GetId();
    }
  }
  Assert(src_edge.value != NoId);
  Assert(dst_edge.value != NoId);

  Assert(not first_src_node.template IsOverlaid<Neighbors>());
  first_src_node.template SetOverlay<Neighbors>();
  auto sc = src_parent_edge.GetClade().value;
  auto dc = dst_parent_edge.GetClade().value;
  Assert(not src_parent_edge.template IsOverlaid<Endpoints>());
  src_parent_edge.template SetOverlay<Endpoints>();
  src_parent_edge.Set(src_parent_node, first_src_node, {sc});
  Assert(not dst_parent_edge.template IsOverlaid<Endpoints>());
  dst_parent_edge.template SetOverlay<Endpoints>();
  dst_parent_edge.Set(dst_parent_node, first_dst_node, {dc});

  first_src_node.SetSingleParent(src_parent_edge);
  Assert(not first_dst_node.template IsOverlaid<Neighbors>());
  first_dst_node.template SetOverlay<Neighbors>();
  first_dst_node.SetSingleParent(dst_parent_edge);

  if (is_sibling_move) {
    auto src_parent_single_parent = src_parent_node.GetParentsCount() > 0
                                        ? src_parent_node.GetSingleParent()
                                        : EdgeId{NoId};
    if (not src_parent_node.template IsOverlaid<Neighbors>()) {
      src_parent_node.template SetOverlay<Neighbors>();
    }
    src_parent_node.ClearConnections();
    if (src_parent_single_parent != EdgeId{NoId}) {
      src_parent_node.SetSingleParent(src_parent_single_parent);
    }
    size_t clade_ctr = 0;
    for (auto e_id : src_sibling_edges) {
      if (e_id != dst_edge) {
        auto e = dag.Get(e_id);
        auto e_child = e.GetChild();
        if (not e.template IsOverlaid<Endpoints>()) {
          e.template SetOverlay<Endpoints>();
        }
        e.Set(src_parent_node, e_child, {clade_ctr});
        src_parent_node.AddEdge({clade_ctr++}, e, true);
        if (not e_child.template IsOverlaid<Neighbors>()) {
          e_child.template SetOverlay<Neighbors>();
        }
        e_child.SetSingleParent(e);
      }
    }
    if (not new_node.IsAppended()) {
      if (not new_node.template IsOverlaid<Neighbors>()) {
        new_node.template SetOverlay<Neighbors>();
      }
    }
    if (not new_edge.IsAppended()) {
      if (not new_edge.template IsOverlaid<Endpoints>()) {
        new_edge.template SetOverlay<Endpoints>();
      }
    }
    new_edge.Set(src_parent_node, new_node, {clade_ctr});
    src_parent_node.AddEdge({clade_ctr}, new_edge, true);
    new_node.SetSingleParent(new_edge);
    src_parent_edge.Set(new_node, first_src_node, {0});
    new_node.AddEdge({0}, src_parent_edge, true);
    first_src_node.SetSingleParent(src_parent_edge);
    dst_parent_edge.Set(new_node, first_dst_node, {1});
    first_dst_node.SetSingleParent(dst_parent_edge);
    new_node.AddEdge({1}, dst_parent_edge, true);
  } else if (has_unifurcation_after_move) {
    auto src_grandparent = src_parent_node.GetSingleParent().GetParent();
    auto dst_parent_single_parent = dst_parent_node.GetParentsCount() > 0
                                        ? dst_parent_node.GetSingleParent()
                                        : EdgeId{NoId};
    auto src_grandparent_single_parent = src_grandparent.GetParentsCount() > 0
                                             ? src_grandparent.GetSingleParent()
                                             : EdgeId{NoId};
    size_t clade_ctr = 0;

    std::vector<EdgeId> src_gp_children;
    for (auto c : src_grandparent.GetChildren()) {
      src_gp_children.push_back(c);
    }
    if (not src_grandparent.template IsOverlaid<Neighbors>()) {
      src_grandparent.template SetOverlay<Neighbors>();
    }
    if (not dst_parent_node.template IsOverlaid<Neighbors>()) {
      dst_parent_node.template SetOverlay<Neighbors>();
    }
    if (not new_edge.template IsOverlaid<Endpoints>()) {
      new_edge.template SetOverlay<Endpoints>();
    }
    if (not new_node.template IsOverlaid<Neighbors>()) {
      new_node.template SetOverlay<Neighbors>();
    }
    std::vector<EdgeId> src_parent_siblings;
    for (auto e_id : src_gp_children) {
      if (e_id != new_edge.GetId() and (e_id != dst_edge)) {
        src_parent_siblings.push_back(e_id);
      }
    }
    src_grandparent.ClearConnections();
    if (src_grandparent_single_parent.value != NoId) {
      src_grandparent.SetSingleParent(src_grandparent_single_parent);
    }
    for (auto e_id : src_parent_siblings) {
      if (e_id != dst_edge) {
        auto e = dag.Get(e_id);
        auto e_child = e.GetChild();
        if (not e.template IsOverlaid<Endpoints>()) {
          e.template SetOverlay<Endpoints>();
        }
        e.Set(src_grandparent, e_child, {clade_ctr});
        src_grandparent.AddEdge({clade_ctr++}, e, true);
        if (not e_child.template IsOverlaid<Neighbors>()) {
          e_child.template SetOverlay<Neighbors>();
        }
        e_child.SetSingleParent(e);
      }
    }
    for (auto e_id : src_sibling_edges) {
      auto e = dag.Get(e_id);
      auto e_child = e.GetChild();
      if (not e.template IsOverlaid<Endpoints>()) {
        e.template SetOverlay<Endpoints>();
      }
      e.Set(src_grandparent, e_child, {clade_ctr});
      src_grandparent.AddEdge({clade_ctr++}, e, true);
      if (not e_child.template IsOverlaid<Neighbors>()) {
        e_child.template SetOverlay<Neighbors>();
      }
      e_child.SetSingleParent(e);
    }
    if (dst_parent_node != src_grandparent) {
      if (not dst_parent_node.template IsOverlaid<Neighbors>()) {
        dst_parent_node.template SetOverlay<Neighbors>();
      }
      dst_parent_node.ClearConnections();
      if (dst_parent_single_parent.value != NoId) {
        auto dst_parent_single_parent_edge = dag.Get(dst_parent_single_parent);
        auto dst_parent_single_parent_clade =
            dst_parent_single_parent_edge.GetClade().value;
        auto p = dst_parent_single_parent_edge.GetParent();
        if (not dst_parent_single_parent_edge.template IsOverlaid<Endpoints>()) {
          dst_parent_single_parent_edge.template SetOverlay<Endpoints>();
        }
        dst_parent_node.SetSingleParent(dst_parent_single_parent_edge);
        dst_parent_single_parent_edge.Set(p, dst_parent_node,
                                          {dst_parent_single_parent_clade});
      }
      clade_ctr = 0;
    }
    for (auto e_id : dst_sibling_edges) {
      if (e_id != new_edge.GetId() and
          std::find(src_parent_siblings.begin(), src_parent_siblings.end(), e_id) ==
              src_parent_siblings.end()) {
        auto e = dag.Get(e_id);
        auto e_child = e.GetChild();
        if (not e.template IsOverlaid<Endpoints>()) {
          e.template SetOverlay<Endpoints>();
        }
        e.Set(dst_parent_node, e_child, {clade_ctr});
        dst_parent_node.AddEdge({clade_ctr++}, e, true);
        if (not e_child.template IsOverlaid<Neighbors>()) {
          e_child.template SetOverlay<Neighbors>();
        }
        e_child.SetSingleParent(e);
      }
    }
    new_node.ClearConnections();
    new_edge.Set(dst_parent_node, new_node, {clade_ctr});
    dst_parent_node.AddEdge({clade_ctr++}, new_edge, true);
    new_node.SetSingleParent(new_edge);
    src_parent_edge.Set(new_node, first_src_node, {0});
    new_node.AddEdge({0}, src_parent_edge, true);
    first_src_node.SetSingleParent(src_parent_edge);
    dst_parent_edge.Set(new_node, first_dst_node, {1});
    new_node.AddEdge({1}, dst_parent_edge, true);
    first_dst_node.SetSingleParent(dst_parent_edge);
  } else {
    auto src_grandparent_edge = src_parent_node.GetParentsCount() > 0
                                    ? src_parent_node.GetSingleParent()
                                    : EdgeId{NoId};
    auto dst_grandparent_edge = dst_parent_node.GetParentsCount() > 0
                                    ? dst_parent_node.GetSingleParent()
                                    : EdgeId{NoId};
    if (not src_parent_node.template IsOverlaid<Neighbors>()) {
      src_parent_node.template SetOverlay<Neighbors>();
    }
    if (not dst_parent_node.template IsOverlaid<Neighbors>()) {
      dst_parent_node.template SetOverlay<Neighbors>();
    }
    src_parent_node.ClearConnections();
    dst_parent_node.ClearConnections();
    if (src_grandparent_edge != EdgeId{NoId}) {
      src_parent_node.SetSingleParent(src_grandparent_edge);
    }
    if (dst_grandparent_edge != EdgeId{NoId}) {
      dst_parent_node.SetSingleParent(dst_grandparent_edge);
    }
    size_t src_parent_clade_ctr = 0;
    for (auto src_sib_edge_id : src_sibling_edges) {
      auto src_sib_edge = dag.Get(src_sib_edge_id);
      auto src_sib_node = src_sib_edge.GetChild();
      if (not src_sib_edge.template IsOverlaid<Endpoints>()) {
        src_sib_edge.template SetOverlay<Endpoints>();
      }
      src_sib_edge.Set(src_parent_node, src_sib_node, {src_parent_clade_ctr});
      src_parent_node.AddEdge({src_parent_clade_ctr++}, src_sib_edge, true);
      if (not src_sib_node.template IsOverlaid<Neighbors>()) {
        src_sib_node.template SetOverlay<Neighbors>();
      }
      src_sib_node.SetSingleParent(src_sib_edge);
    }
    size_t dst_parent_clade_ctr = 0;
    for (auto dst_sib_edge_id : dst_sibling_edges) {
      auto dst_sib_edge = dag.Get(dst_sib_edge_id);
      auto dst_sib_node = dst_sib_edge.GetChild();
      if (not dst_sib_edge.template IsOverlaid<Endpoints>()) {
        dst_sib_edge.template SetOverlay<Endpoints>();
      }
      dst_sib_edge.Set(dst_parent_node, dst_sib_node, {dst_parent_clade_ctr});
      dst_parent_node.AddEdge({dst_parent_clade_ctr++}, dst_sib_edge, true);
      if (not dst_sib_node.template IsOverlaid<Neighbors>()) {
        dst_sib_node.template SetOverlay<Neighbors>();
      }
      dst_sib_node.SetSingleParent(dst_sib_edge);
    }
    if (not new_edge.IsAppended()) {
      if (not new_edge.template IsOverlaid<Endpoints>()) {
        new_edge.template SetOverlay<Endpoints>();
      }
    }
    if (not new_node.IsAppended()) {
      if (not new_node.template IsOverlaid<Neighbors>()) {
        new_node.template SetOverlay<Neighbors>();
      }
      new_node.ClearConnections();
    }
    src_parent_edge.Set(new_node, first_src_node, {0});
    new_node.AddEdge({0}, src_parent_edge, true);
    first_src_node.SetSingleParent(src_parent_edge);
    dst_parent_edge.Set(new_node, first_dst_node, {1});
    new_node.AddEdge({1}, dst_parent_edge, true);
    first_dst_node.SetSingleParent(dst_parent_edge);
    new_edge.Set(dst_parent_node, new_node, {dst_parent_clade_ctr});
    dst_parent_node.AddEdge({dst_parent_clade_ctr}, new_edge, true);
    new_node.SetSingleParent(new_edge);
  }
  Assert(new_node.GetParentsCount() == 1);
  Assert(dag.GetNodesCount() ==
         (has_unifurcation_after_move ? old_num_nodes : old_num_nodes + 1));
  Assert(dag.GetEdgesCount() ==
         (has_unifurcation_after_move ? old_num_edges : old_num_edges + 1));
  return {new_node, has_unifurcation_after_move};
}

}  // namespace

template <typename DAG, typename Backend, typename CRTP, typename Tag>
std::pair<NodeId, bool> FeatureMutableView<HypotheticalTree<DAG, Backend>, CRTP,
                                           Tag>::ApplyMove(NodeId lca, NodeId src,
                                                           NodeId dst) const {
  auto& dag = static_cast<const CRTP&>(*this);
  return ApplyMoveImpl(dag, lca, src, dst);
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
bool FeatureMutableView<HypotheticalTree<DAG, Backend>, CRTP, Tag>::
    InitHypotheticalTree(const Profitable_Moves& move,
                         const std::vector<Node_With_Major_Allele_Set_Change>&
                             nodes_with_major_allele_set_change) const {
  auto& self = GetFeatureStorage(this).get();
  Assert(not self.data_);
  auto& dag = static_cast<const CRTP&>(*this);

  auto [new_node, has_unifurcation_after_move] =
      dag.ApplyMove(dag.GetNodeFromMAT(move.LCA), dag.GetNodeFromMAT(move.src),
                    dag.GetNodeFromMAT(move.dst));

  if (new_node.value == NoId) {
    return false;
  }

  self.data_ = std::make_unique<typename HypotheticalTree<DAG, Backend>::Data>(
      dag, move, new_node, has_unifurcation_after_move,
      nodes_with_major_allele_set_change);

  if (dag.GetMoveLCA().IsUA()) {
    return true;
  }
  NodeId lca = dag.GetMoveLCA();
  do {
    self.data_->lca_ancestors_.insert(lca);
    lca = dag.Get(lca).GetSingleParent().GetParent();
  } while (not dag.Get(lca).IsUA());

  auto mark_changed = [&dag](NodeId current_node) {
    while (not(current_node == dag.GetMoveLCA().GetId() or
               dag.Get(current_node).HasChangedTopology())) {
      if (not dag.Get(current_node).IsAppended()) {
        dag.Get(current_node).template SetOverlay<HypotheticalNode>();
      }
      dag.Get(current_node).SetHasChangedTopology();
      if (dag.Get(current_node).IsUA()) {
        break;
      }
      current_node = dag.Get(current_node).GetSingleParent().GetParent();
    }
  };

  mark_changed(dag.GetOldSourceParent());
  mark_changed(dag.GetMoveTarget().GetOld().GetSingleParent().GetParent());
  mark_changed(dag.GetNodeFromMAT(move.dst));
  mark_changed(dag.GetNodeFromMAT(move.src));
  mark_changed(dag.GetNodeFromMAT(move.dst->parent));
  mark_changed(dag.GetNodeFromMAT(move.src->parent));
  // for the case of a unifurcation after move, we need to check both pre-move
  // grandparents
  if (move.dst->parent != nullptr) {
    if (move.dst->parent->parent != nullptr) {
      mark_changed(dag.GetNodeFromMAT(move.dst->parent->parent));
    }
  }
  if (move.src->parent != nullptr) {
    if (move.src->parent->parent != nullptr) {
      mark_changed(dag.GetNodeFromMAT(move.src->parent->parent));
    }
  }
  mark_changed(dag.GetMoveNew());

  return true;
}

template <typename DAG, typename Backend, typename CRTP, typename Tag>
bool FeatureMutableView<HypotheticalTree<DAG, Backend>, CRTP,
                        Tag>::InitHypotheticalTree(NodeId src, NodeId dst,
                                                   NodeId lca) const {
  auto& self = GetFeatureStorage(this).get();
  Assert(not self.data_);
  auto& dag = static_cast<const CRTP&>(*this);

  auto [new_node, has_unifurcation_after_move] = dag.ApplyMove(lca, src, dst);

  if (new_node.value == NoId) {
    return false;
  }

  self.data_ = std::make_unique<typename HypotheticalTree<DAG, Backend>::Data>(
      dag, src, dst, lca, new_node, has_unifurcation_after_move);

  if (dag.GetMoveLCA().IsUA()) {
    return true;
  }
  NodeId lca_id = dag.GetMoveLCA();
  do {
    self.data_->lca_ancestors_.insert(lca_id);
    lca_id = dag.Get(lca_id).GetSingleParent().GetParent();
  } while (not dag.Get(lca_id).IsUA());

  auto mark_changed = [&dag](NodeId current_node) {
    while (not(current_node == dag.GetMoveLCA().GetId() or
               dag.Get(current_node).HasChangedTopology())) {
      if (not dag.Get(current_node).IsAppended()) {
        dag.Get(current_node).template SetOverlay<HypotheticalNode>();
      }
      dag.Get(current_node).SetHasChangedTopology();
      if (dag.Get(current_node).IsUA()) {
        break;
      }
      current_node = dag.Get(current_node).GetSingleParent().GetParent();
    }
  };

  mark_changed(dag.GetOldSourceParent());
  mark_changed(dag.GetMoveTarget().GetOld().GetSingleParent().GetParent());
  mark_changed(dst);
  mark_changed(src);
  mark_changed(dag.Get(dst).GetOld().GetSingleParent().GetParent());
  mark_changed(dag.Get(src).GetOld().GetSingleParent().GetParent());
  // for the case of a unifurcation after move, check pre-move grandparents
  auto dst_parent = dag.Get(dst).GetOld().GetSingleParent().GetParent();
  if (not dst_parent.IsUA()) {
    mark_changed(dst_parent.GetSingleParent().GetParent());
  }
  auto src_parent = dag.Get(src).GetOld().GetSingleParent().GetParent();
  if (not src_parent.IsUA()) {
    mark_changed(src_parent.GetSingleParent().GetParent());
  }
  mark_changed(dag.GetMoveNew());

  return true;
}

// HypotheticalTree::Data constructors

template <typename DAG, typename Backend>
template <typename DAGView>
HypotheticalTree<DAG, Backend>::Data::Data(
    const DAGView& dag, const Profitable_Moves& move, NodeId new_node,
    bool has_unifurcation_after_move,
    const std::vector<Node_With_Major_Allele_Set_Change>&
        nodes_with_major_allele_set_change)
    : new_node_{new_node},
      has_unifurcation_after_move_{has_unifurcation_after_move},
      move_{move} {
  backend_.Initialize(dag, move, nodes_with_major_allele_set_change);
}

template <typename DAG, typename Backend>
template <typename DAGView>
HypotheticalTree<DAG, Backend>::Data::Data(const DAGView& dag, NodeId src, NodeId dst,
                                           NodeId lca, NodeId new_node,
                                           bool has_unifurcation_after_move)
    : new_node_{new_node},
      has_unifurcation_after_move_{has_unifurcation_after_move},
      move_{} {
  // Detect MatOptimizeScoringBackend by checking if it has
  // Initialize(dag, Profitable_Moves, vector<Node_With_Major_Allele_Set_Change>).
  // All other backends (ML, ParsimonyOnly) use NodeId-based Initialize.
  if constexpr (requires(Backend b, DAGView d, Profitable_Moves m,
                         std::vector<Node_With_Major_Allele_Set_Change> c) {
                  b.Initialize(d, m, c);
                }) {
    // matOptimize backend: convert NodeIds to MAT pointers
    move_.src = dag.Get(src).GetMATNode();
    move_.dst = dag.Get(dst).GetMATNode();
    move_.LCA = dag.Get(lca).GetMATNode();
    move_.score_change = 0;
    std::vector<Node_With_Major_Allele_Set_Change> empty_changes;
    backend_.Initialize(dag, move_, empty_changes);
  } else {
    backend_.Initialize(dag, src, dst, lca);
  }
}

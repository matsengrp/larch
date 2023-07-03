template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsMATRoot() const {
  auto& node = static_cast<const CRTP&>(*this);
  return node.GetMATNode()->parent == nullptr;
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsMoveSource() const {
  auto& node = static_cast<const CRTP&>(*this);
  auto move_src = node.GetDAG().GetMoveSources();
  //return node.GetDAG().GetMoveSource().GetId() == node.GetId();
  return (std::find(move_src.begin(), move_src.end(), node.GetId()) != move_src.end());
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsMoveTarget() const {
  auto& node = static_cast<const CRTP&>(*this);
  auto move_dst = node.GetDAG().GetMoveTargets();
  //return node.GetDAG().GetMoveTarget().GetId() == node.GetId();
  return (std::find(move_dst.begin(), move_dst.end(), node.GetId()) != move_dst.end());
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
  return GetFeatureStorage(this).has_changed_topology_;
}

template <typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalNode, CRTP, Tag>::GetOld() const {
  auto& node = static_cast<const CRTP&>(*this);
  return node.GetDAG().GetOriginal().Get(node.GetId());
}

template <typename CRTP, typename Tag>
const ContiguousSet<MutationPosition>&
FeatureConstView<HypotheticalNode, CRTP, Tag>::GetChangedBaseSites() const {
  return GetFeatureStorage(this).changed_base_sites_;
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
  if (node.GetOld().IsLeaf() or node.IsMoveTarget()) {
    // if it's a leaf node, then the fitch sets don't change.
    // if it's the target node, then fitch sets can't have changed,
    // but the fitch sets recorded in tree_'s changed fitch set map
    // relative to this node are meant for this node's new parent!
    return {node.GetMATNode()->mutations, std::nullopt};
  } else if (node.IsMoveNew()) {
    // Then fitch set changes are relative to the target node
    auto result = dag.GetChangedFitchSetMap().find(dag.GetMoveTarget().GetMATNode());
    return {dag.GetMoveTarget().GetMATNode()->mutations,
            result == dag.GetChangedFitchSetMap().end()
                ? std::nullopt
                : std::make_optional(result->second.Copy())};
  } else {
    auto result = dag.GetChangedFitchSetMap().find(node.GetMATNode());
    return {node.GetMATNode()->mutations,
            result == dag.GetChangedFitchSetMap().end()
                ? std::nullopt
                : std::make_optional(result->second.Copy())};
  }
}

template <typename CRTP, typename Tag>
FitchSet FeatureConstView<HypotheticalNode, CRTP, Tag>::GetFitchSetAtSite(
    MutationPosition site) const {
  auto node = static_cast<const CRTP&>(*this).Const();
  auto dag = node.GetDAG();
  Assert(site.value <= dag.GetReferenceSequence().size());
  auto [old_fitch_sets, changes] =
      GetFitchSetParts();  // TODO: modify to also return boundary alleles

  auto build_new_fitch_set =
      [&node](nuc_one_hot old_major_alleles, nuc_one_hot old_boundary_alleles,
              nuc_one_hot rem_set, nuc_one_hot add_set) -> FitchSet {
    Assert(old_major_alleles <= 16);
    Assert(old_boundary_alleles <= 16);
    Assert(rem_set <= 16);
    Assert(add_set <= 16);
    if (old_major_alleles & add_set) {
      auto result = old_major_alleles & add_set;
      Assert(0 < result);
      Assert(result <= 16);
      return FitchSet(result);
    } else {
      auto result = (old_major_alleles & ~rem_set) | (old_boundary_alleles & add_set);
      Assert(result <= 16);
      if (result > 0) {
        return FitchSet(result);
      } else {
        result = old_major_alleles | (old_boundary_alleles & ~rem_set);
        return FitchSet(result);
      }
    }
  };

  if (old_fitch_sets.find(static_cast<int>(site.value)) ==
      old_fitch_sets.mutations.end()) {
    // if no fitch set is recorded on the corresponding MAT node, we can use
    // a singleton set containing the base at this site in the parent compact genome
    // In the special case that the node is the new move, its Fitch sets (and Fitch set
    // changes) should be calculated relative to the old parent of GetMoveTarget, and so
    // we will retrieve the base of that parent for the Fitch set in this case.

    nuc_one_hot boundary_allele =
        0;  // from the 3 conditions for empty fitch set (see Cheng's comment)
    auto old_parent_base =
        node.IsMoveNew()
            ? base_to_singleton(dag.GetMoveTarget()
                                    .GetOld()
                                    .GetSingleParent()
                                    .GetParent()
                                    .GetCompactGenome()
                                    .GetBase(site, dag.GetReferenceSequence()))
            : base_to_singleton(node.GetOld()
                                    .GetSingleParent()
                                    .GetParent()
                                    .GetCompactGenome()
                                    .GetBase(site, dag.GetReferenceSequence()));
    if (changes.has_value() and changes.value().Contains(site)) {
      return build_new_fitch_set(old_parent_base, boundary_allele,
                                 changes.value().at(site).get_decremented(),
                                 changes.value().at(site).get_incremented());
    } else {
      return FitchSet(old_parent_base);
    }
  } else if (changes.has_value() and changes.value().Contains(site)) {
    return FitchSet(build_new_fitch_set(
        old_fitch_sets.find(static_cast<int>(site.value))->get_all_major_allele(),
        old_fitch_sets.find(static_cast<int>(site.value))->get_boundary1_one_hot(),
        changes.value().at(site).get_decremented(),
        changes.value().at(site).get_incremented()));
  } else {
    return FitchSet(
        old_fitch_sets.find(static_cast<int>(site.value))->get_all_major_allele());
  }
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
  Assert(node.HaveMATNode());
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
  CompactGenome result = old_cg.Copy();
  result.ApplyChanges(cg_changes);
  return result;
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsNonrootAnchorNode() const {
  auto node = static_cast<const CRTP&>(*this).Const();
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
  GetFeatureStorage(this).has_changed_topology_ = true;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<HypotheticalNode, CRTP, Tag>::PreorderComputeCompactGenome(
    std::vector<NodeId>& result_nodes, std::vector<EdgeId>& result_edges) const {
  auto& node = static_cast<const CRTP&>(*this);
  if (not node.IsUA() and not node.IsMoveNew()) {
    node.template SetOverlay<Deduplicate<CompactGenome>>();
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

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetMoveLCA() const {
  auto& self = GetFeatureStorage(this);
  Assert(self.data_);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetNodeFromMAT(self.data_->move_.LCA);
  //return dag.Get(dag.GetUncondensedNodeFromMAT(self.data_->move_.LCA)[0]);
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetMoveSource() const {
  auto& self = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetNodeFromMAT(self.data_->move_.src);
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetMoveSources() const {
  auto& self = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetUncondensedNodeFromMAT(self.data_->move_.src);
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetMoveTarget() const {
  auto& self = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetNodeFromMAT(self.data_->move_.dst);
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetMoveTargets() const {
  auto& self = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetUncondensedNodeFromMAT(self.data_->move_.dst);
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetMoveNew() const {
  auto& self = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.Get(self.data_->new_node_);
}

template <typename DAG, typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::HasUnifurcationAfterMove()
    const {
  auto& self = GetFeatureStorage(this);
  return self.data_->has_unifurcation_after_move_;
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetOldSourceParent() const {
  auto& dag = static_cast<const CRTP&>(*this);
  if (dag.HasUnifurcationAfterMove()) {
    return dag.GetMoveSource()
        .GetOld()
        .GetSingleParent()
        .GetParent()
        .GetSingleParent()
        .GetParent();
  } else {
    return dag.GetMoveSource().GetOld().GetSingleParent().GetParent();
  }
}

template <typename Node>
static inline bool IsChanged(Node node) {
  return node.GetCompactGenome() != node.GetOld().GetCompactGenome();
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetOldestChangedNode() const {
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
      auto srcs = dag.GetMoveSources();
      if (std::find(srcs.begin(), srcs.end(), node.GetId()) != srcs.end()) {
      //if (node.GetId() != dag.GetMoveSource().GetId()) {
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

template <typename DAG, typename CRTP, typename Tag>
std::pair<std::vector<NodeId>, std::vector<EdgeId>>
FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetFragment() const {
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

  return {result_nodes, result_edges};
//  auto collapsed = dag.CollapseEmptyFragmentEdges(result_nodes, result_edges);
//  return {collapsed.first, collapsed.second};
}

template <typename DAG, typename CRTP, typename Tag>
std::pair<std::vector<NodeId>, std::vector<EdgeId>>
FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::CollapseEmptyFragmentEdges(
    std::vector<NodeId> fragment_nodes, std::vector<EdgeId> fragment_edges) const {
  auto& dag = static_cast<const CRTP&>(*this);

  // fragment_nodes is computed in a preorder traversal
  auto fragment_root = fragment_nodes.front();

  // keep track of edges/nodes that are collapsible
  std::unordered_map<NodeId, bool> is_parent_of_collapsible_edge;
  std::unordered_map<NodeId, bool> is_child_of_collapsible_edge;
  std::unordered_map<EdgeId, bool> is_collapsible_edge;
  for (auto edge_id : fragment_edges) {
    auto edge = dag.Get(edge_id);
    auto parent = edge.GetParent();
    auto child = edge.GetChild();
    if (not(parent.Const().GetCompactGenome() != child.Const().GetCompactGenome() or
            child.IsLeaf() or parent.IsUA())) {
      is_parent_of_collapsible_edge.insert({parent, true});
      is_child_of_collapsible_edge.insert({child, true});
      is_collapsible_edge.insert({edge, true});
    }
  }

  std::vector<NodeId> current_nodes;
  std::vector<EdgeId> current_edges;

  for (auto node_id : fragment_nodes) {
    if (not is_child_of_collapsible_edge[node_id]) {
      current_nodes.push_back(node_id);
      if (is_parent_of_collapsible_edge[node_id]) {
        auto parent = dag.Get(node_id);

        std::vector<EdgeId> current_children;
        if (parent.IsMoveNew()) {
          current_children.push_back(dag.GetMoveSource().GetSingleParent());
          current_children.push_back(dag.GetMoveTarget().GetSingleParent());
        } else {
          for (auto c : parent.GetChildren() | Transform::GetId()) {
            current_children.push_back(c);
          }
        }
        auto grandparent_edge = parent.GetSingleParent();
        size_t current_clade = 0;

        parent.template SetOverlay<Neighbors>();
        parent.ClearConnections();

        bool collapsed_all_children = false;
        while (not collapsed_all_children) {
          collapsed_all_children = true;
          std::vector<EdgeId> still_to_collapse;

          // try to add all (non-collapsible) edges below current node, as well as the
          // children of any collapsible edges
          for (auto edge_id : current_children) {
            auto edge = dag.Get(edge_id);
            auto child = edge.GetChild();
            if (is_collapsible_edge[edge]) {
              if (child.IsNonrootAnchorNode()) {
                child.template SetOverlay<HypotheticalNode>();
              } else if (node_id == fragment_root) {
                // that is, if parent is root anchor node...
                current_nodes.insert(current_nodes.begin(),
                                     grandparent_edge.GetParent().GetId());
              }
              for (auto child_edge : child.GetChildren()) {
                if (is_collapsible_edge[child_edge]) {
                  still_to_collapse.push_back(child_edge);
                } else {
                  auto grandchild_node = child_edge.GetChild();
                  child_edge.template SetOverlay<Endpoints>();
                  child_edge.Set(parent, grandchild_node, {current_clade});
                  parent.AddEdge({current_clade++}, child_edge, true);
                }
              }
            } else {
              edge.template SetOverlay<Endpoints>();
              edge.Set(parent, child, {current_clade});
              parent.AddEdge({current_clade++}, edge, true);
            }
          }
          // check if any of the grandchildren/edges we're trying to add are in turn
          // collapsible
          if (still_to_collapse.size() > 0) {
            current_children.clear();
            collapsed_all_children = false;
            current_children.insert(current_children.begin(), still_to_collapse.begin(),
                                    still_to_collapse.end());
          }
        }
      }
    }
  }
  for (auto node_id : current_nodes) {
    auto node = dag.Get(node_id);
    for (auto child : node.GetChildren()) {
      current_edges.push_back(child);
    }
  }
  return {current_nodes, current_edges};
}

template <typename DAG, typename CRTP, typename Tag>
const ContiguousMap<MATNodePtr, ContiguousMap<MutationPosition, Mutation_Count_Change>>&
FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetChangedFitchSetMap() const {
  auto& self = GetFeatureStorage(this);
  Assert(self.data_);
  return self.data_->changed_fitch_set_map_;
}

template <typename DAG, typename CRTP, typename Tag>
const ContiguousSet<NodeId>&
FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetLCAAncestors() const {
  auto& self = GetFeatureStorage(this);
  Assert(self.data_);
  return self.data_->lca_ancestors_;
}

namespace {

template <typename DAG>
std::pair<NodeId, bool> ApplyMoveImpl(DAG dag, std::vector<NodeId>& src, std::vector<NodeId>& dst) {
  for (auto node: dag.GetNodes()) {
    if (not node.IsUA() and not node.IsMATRoot()) {
     Assert(node.GetParentsCount() == 1);
    }
  }
  Assert(dag.IsTree());

  auto first_src_node = dag.Get(src[0]);
  auto first_dst_node = dag.Get(dst[0]);

  if (first_src_node.IsTreeRoot() or first_src_node.GetId() == first_dst_node.GetId() or first_dst_node.IsTreeRoot()) {
    // no-op
    return {};
  }

  auto src_parent_node = first_src_node.GetSingleParent().GetParent();
  auto dst_parent_node = first_dst_node.GetSingleParent().GetParent();
  const bool is_sibling_move = src_parent_node.GetId() == dst_parent_node.GetId();

  const bool has_unifurcation_after_move = src_parent_node.GetCladesCount() == src.size() + 1;
  if ((is_sibling_move and has_unifurcation_after_move) or src_parent_node.IsTreeRoot()) {
    // no-op
    return {};
  }

  auto new_node = has_unifurcation_after_move ? src_parent_node : dag.AppendNode();
  auto new_edge = has_unifurcation_after_move ? src_parent_node.GetSingleParent() : dag.AppendEdge();

  std::vector<EdgeId> src_edges;
  std::vector<EdgeId> src_sibling_edges;
  for (auto e: src_parent_node.GetChildren()) {
    if (std::find(src.begin(), src.end(), e.GetChild().GetId()) == src.end()) {
      src_sibling_edges.push_back(e.GetId());
    } else {
      src_edges.push_back(e.GetId());
    }
  }
  std::vector<EdgeId> dst_edges;
  std::vector<EdgeId> dst_sibling_edges;
  for (auto e: dst_parent_node.GetChildren()) {
    if (std::find(dst.begin(), dst.end(), e.GetChild().GetId()) == dst.end()) {
      dst_sibling_edges.push_back(e.GetId());
    } else {
      dst_edges.push_back(e.GetId());
    }
  }
  Assert(src.size() == src_edges.size());
  Assert(dst.size() == dst_edges.size());
  std::cout << "applying move " << first_src_node.GetId().value << " -> " << first_dst_node.GetId().value << "\n" << std::flush;

  if (is_sibling_move) {
    std::cout << "type 3\n" << std::flush;
    auto src_parent_single_parent = src_parent_node.GetParentsCount() > 0 ? src_parent_node.GetSingleParent() : EdgeId{NoId};
    src_parent_node.template SetOverlay<Neighbors>();
    new_node.template SetOverlay<Neighbors>();
    new_edge.template SetOverlay<Endpoints>();
    src_parent_node.ClearConnections();
    src_parent_node.SetSingleParent(src_parent_single_parent);
    new_node.ClearConnections();
    size_t clade_ctr = 0;
    for (auto e_id: src_sibling_edges) {
      if (std::find(dst_edges.begin(), dst_edges.end(), e_id) == dst_edges.end()) {
        auto e = dag.Get(e_id);
        e.template SetOverlay<Endpoints>();
        e.Set(src_parent_node, e.GetChild(), {clade_ctr});
        src_parent_node.AddEdge({clade_ctr++}, e, true);
      }
    }
    new_edge.Set(src_parent_node, new_node, {clade_ctr});
    src_parent_node.AddEdge({clade_ctr}, new_edge, true);
    clade_ctr = 0;
    for (auto e_id: src_edges) {
      auto e = dag.Get(e_id);
      e.template SetOverlay<Endpoints>();
      e.Set(new_node, e.GetChild(), {clade_ctr});
      new_node.AddEdge({clade_ctr++}, e, true);
    }
    for (auto e_id: dst_edges) {
      auto e = dag.Get(e_id);
      e.template SetOverlay<Endpoints>();
      e.Set(new_node, e.GetChild(), {clade_ctr});
      new_node.AddEdge({clade_ctr++}, e, true);
    }
    for (auto c: src_parent_node.GetChildren()) {
      std::cout << "src parent: " << src_parent_node.GetId().value << " - " << c.GetChild().GetId().value << "\n";
    }
    for (auto c: new_node.GetChildren()) {
      std::cout << "new_node  : " << new_node.GetId().value << " - " << c.GetChild().GetId().value << "\n";
    }
  } else if (has_unifurcation_after_move) {
    std::cout << "type 1\n" << std::flush;
    auto src_grandparent = src_parent_node.GetSingleParent().GetParent();
    auto src_sib_edge = dag.Get(src_sibling_edges.at(0));
    auto src_sib_node = src_sib_edge.GetChild();
    auto dst_parent_single_parent = dst_parent_node.GetParentsCount() > 0 ? dst_parent_node.GetSingleParent() : EdgeId{NoId};
    auto src_grandparent_single_parent = src_grandparent.GetParentsCount() > 0 ? src_grandparent.GetSingleParent() : EdgeId{NoId};
    size_t clade_ctr = 0;

    src_grandparent.template SetOverlay<Neighbors>();
    dst_parent_node.template SetOverlay<Neighbors>();
    new_node.template SetOverlay<Neighbors>();
    new_edge.template SetOverlay<Endpoints>();

    auto src_grandparent_edge = src_parent_node.GetSingleParent();
    src_parent_node.ClearConnections();
    src_parent_node.SetSingleParent(new_edge);
    std::vector<EdgeId> src_parent_siblings;
    for (auto e: src_grandparent.GetChildren()) {
      if (e.GetId() != src_grandparent_edge.GetId() and (std::find(dst_edges.begin(), dst_edges.end(), e.GetId()) == dst_edges.end())) {
        src_parent_siblings.push_back(e.GetId());
      }
    }
    src_grandparent.ClearConnections();
    if (src_grandparent_single_parent.value != NoId) {
      src_grandparent.SetSingleParent(src_grandparent_single_parent);
    }
    for (auto e_id: src_parent_siblings) {
      if (std::find(dst_edges.begin(), dst_edges.end(), e_id) == dst_edges.end()) {
        auto e = dag.Get(e_id);
        e.template SetOverlay<Endpoints>();
        e.Set(src_grandparent, e.GetChild(), {clade_ctr});
        src_grandparent.AddEdge({clade_ctr++}, e, true);
std::cout << "src gp1: " << src_grandparent.GetId().value << " -> " << e.GetChild().GetId().value << "\n";
      }
    }
    for (auto e_id: src_sibling_edges) {
      auto e = dag.Get(e_id);
      e.template SetOverlay<Endpoints>();
      e.Set(src_grandparent, e.GetChild(), {clade_ctr});
      src_grandparent.AddEdge({clade_ctr++}, e, true);
std::cout << "src gp2: " << src_grandparent.GetId().value << " -> " << e.GetChild().GetId().value << "\n";
    }
    if (dst_parent_node != src_grandparent) {
      dst_parent_node.ClearConnections();
      if (dst_parent_single_parent.value != NoId) {
        dst_parent_node.SetSingleParent(dst_parent_single_parent);
      }
      clade_ctr = 0;
    }
    for (auto e_id: dst_sibling_edges) {
      if (e_id != new_edge.GetId()) {
        auto e = dag.Get(e_id);
        e.template SetOverlay<Endpoints>();
        e.Set(dst_parent_node, e.GetChild(), {clade_ctr});
        dst_parent_node.AddEdge({clade_ctr++}, e, true);
std::cout << "dst  p1: " << dst_parent_node.GetId().value << " -> " << e.GetChild().GetId().value << "\n";
      }
    }
    new_edge.template SetOverlay<Endpoints>();
    new_edge.Set(dst_parent_node, new_node, {clade_ctr});
    dst_parent_node.AddEdge({clade_ctr++}, new_edge, true);
    clade_ctr = 0;
    for (auto e_id: src_edges) {
      auto e = dag.Get(e_id);
      e.template SetOverlay<Endpoints>();
      e.Set(new_node, e.GetChild(), {clade_ctr});
      new_node.AddEdge({clade_ctr++}, e, true);
std::cout << "newnode1: " << new_node.GetId().value << " -> " << e.GetChild().GetId().value << "\n";
    }
    for (auto e_id: dst_edges) {
      auto e = dag.Get(e_id);
      e.template SetOverlay<Endpoints>();
      e.Set(new_node, e.GetChild(), {clade_ctr});
      new_node.AddEdge({clade_ctr++}, e, true);
std::cout << "newnode2: " << new_node.GetId().value << " -> " << e.GetChild().GetId().value << "\n";
    }

for (auto c: src_grandparent.GetChildren()) {
  std::cout << "src gp: " << src_grandparent.GetId().value << " - " << c.GetChild().GetId().value << "\n";
}
for (auto c: dst_parent_node.GetChildren()) {
  std::cout << "dst p : " << dst_parent_node.GetId().value << " - " << c.GetChild().GetId().value << "\n";
}
for (auto c: src_parent_node.GetChildren()) {
  std::cout << "src p : " << src_parent_node.GetId().value << " - " << c.GetChild().GetId().value << "\n";
}
  } else {
    std::cout << "type 2\n" << std::flush;
    auto src_parent_clades_size = src_parent_node.GetCladesCount();
    auto dst_parent_clades_size = dst_parent_node.GetCladesCount();
    src_parent_node.template SetOverlay<Neighbors>();
    dst_parent_node.template SetOverlay<Neighbors>();
    new_node.template SetOverlay<Neighbors>();
    auto src_grandparent_edge = src_parent_node.GetSingleParent();
    auto dst_grandparent_edge = dst_parent_node.GetSingleParent();
    src_parent_node.ClearConnections();
    dst_parent_node.ClearConnections();
    src_parent_node.SetSingleParent(src_grandparent_edge);
    dst_parent_node.SetSingleParent(dst_grandparent_edge);
    new_node.ClearConnections();
    size_t src_parent_clade_ctr = 0;
    for (auto src_sib_edge_id: src_sibling_edges) {
      auto src_sib_edge = dag.Get(src_sib_edge_id);
      auto src_sib_node = src_sib_edge.GetChild();
      src_sib_edge.template SetOverlay<Endpoints>();
      src_sib_edge.Set(src_parent_node, src_sib_node, {src_parent_clade_ctr});
      src_parent_node.AddEdge({src_parent_clade_ctr++}, src_sib_edge, true);
    }
    size_t dst_parent_clade_ctr = 0;
    for (auto dst_sib_edge_id: dst_sibling_edges) {
      auto dst_sib_edge = dag.Get(dst_sib_edge_id);
      auto dst_sib_node = dst_sib_edge.GetChild();
      dst_sib_edge.template SetOverlay<Endpoints>();
      dst_sib_edge.Set(dst_parent_node, dst_sib_node, {dst_parent_clade_ctr});
      dst_parent_node.AddEdge({dst_parent_clade_ctr++}, dst_sib_edge, true);
    }
    size_t new_node_clade_ctr = 0;
    for (auto src_edge_id: src_edges) {
      auto src_edge = dag.Get(src_edge_id);
      auto src_node = src_edge.GetChild();
      src_edge.template SetOverlay<Endpoints>();
      src_edge.Set(new_node, src_node, {new_node_clade_ctr});
      new_node.AddEdge({new_node_clade_ctr++}, src_edge, true);
    }
    for (auto dst_edge_id: dst_edges) {
      auto dst_edge = dag.Get(dst_edge_id);
      auto dst_node = dst_edge.GetChild();
      dst_edge.template SetOverlay<Endpoints>();
      dst_edge.Set(new_node, dst_node, {new_node_clade_ctr});
      new_node.AddEdge({new_node_clade_ctr++}, dst_edge, true);
    }
    new_edge.template SetOverlay<Endpoints>();
    new_edge.Set(dst_parent_node, new_node, {dst_parent_clade_ctr});
    dst_parent_node.AddEdge({dst_parent_clade_ctr}, new_edge, true);
    new_node.SetSingleParent(new_edge);
  }
  return {new_node, has_unifurcation_after_move};
}

}  // namespace

template <typename DAG, typename CRTP, typename Tag>
std::pair<NodeId, bool> FeatureMutableView<HypotheticalTree<DAG>, CRTP, Tag>::ApplyMove(
    NodeId src, NodeId dst) const {
  auto& dag = static_cast<const CRTP&>(*this);
  std::vector<NodeId> srcs{src};
  std::vector<NodeId> dsts{dst};
  return ApplyMoveImpl(dag, srcs, dsts);
}

template <typename DAG, typename CRTP, typename Tag>
std::pair<NodeId, bool> FeatureMutableView<HypotheticalTree<DAG>, CRTP, Tag>::ApplyMove(
    std::vector<NodeId> src, std::vector<NodeId> dst) const {
  auto& dag = static_cast<const CRTP&>(*this);
  return ApplyMoveImpl(dag, src, dst);
}

template <typename DAG, typename CRTP, typename Tag>
bool FeatureMutableView<HypotheticalTree<DAG>, CRTP, Tag>::InitHypotheticalTree(
    const Profitable_Moves& move, const std::vector<Node_With_Major_Allele_Set_Change>&
                                      nodes_with_major_allele_set_change) {
  auto& self = GetFeatureStorage(this);
  Assert(not self.data_);
  auto& dag = static_cast<const CRTP&>(*this);

std::cout << "about to apply move\n" << std::flush;
  auto [new_node, has_unifurcation_after_move] =
      dag.ApplyMove(dag.GetUncondensedNodeFromMAT(move.src), dag.GetUncondensedNodeFromMAT(move.dst));

std::cout << "applied move\n" << std::flush;
  if (new_node.value == NoId) {
    return false;
  }

  self.data_ = std::make_unique<typename HypotheticalTree<DAG>::Data>(
      typename HypotheticalTree<DAG>::Data{move, new_node, has_unifurcation_after_move,
                                           nodes_with_major_allele_set_change});

std::cout << "created data\n" << std::flush;
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
      dag.Get(current_node).template SetOverlay<HypotheticalNode>();
      dag.Get(current_node).SetHasChangedTopology();
      if (dag.Get(current_node).IsUA()) {
        break;
      }
      current_node = dag.Get(current_node).GetSingleParent().GetParent();
    }
  };
  mark_changed(dag.GetOldSourceParent());
  mark_changed(dag.GetMoveNew());

  return true;
}

template <typename DAG>
HypotheticalTree<DAG>::Data::Data(const Profitable_Moves& move, NodeId new_node,
                                  bool has_unifurcation_after_move,
                                  const std::vector<Node_With_Major_Allele_Set_Change>&
                                      nodes_with_major_allele_set_change)
    : move_{move},
      new_node_{new_node},
      has_unifurcation_after_move_{has_unifurcation_after_move} {
  for (auto& node_with_allele_set_change : nodes_with_major_allele_set_change) {
    if (not node_with_allele_set_change.node->is_leaf()) {
      Assert(node_with_allele_set_change.node != nullptr);
      ContiguousMap<MutationPosition, Mutation_Count_Change> node_map;
      for (auto& mutation_count_change :
           node_with_allele_set_change.major_allele_set_change) {
        if (mutation_count_change.get_position() >= 2147483647) {
          continue;
        }
        MutationPosition pos = {
            static_cast<size_t>(mutation_count_change.get_position())};
        node_map.insert({pos, mutation_count_change});
      }
      changed_fitch_set_map_.insert(
          {node_with_allele_set_change.node, std::move(node_map)});
    }
  }
}

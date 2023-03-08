template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsMATRoot() const {
  auto& node = static_cast<const CRTP&>(*this);
  return node.GetMATNode().parent == nullptr;
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsMoveSource() const {
  auto& node = static_cast<const CRTP&>(*this);
  return node.GetDAG().GetMoveSource().GetId() == node.GetId();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsMoveTarget() const {
  auto& node = static_cast<const CRTP&>(*this);
  return node.GetDAG().GetMoveTarget().GetId() == node.GetId();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsMoveNew() const {
  auto& node = static_cast<const CRTP&>(*this);
  return node.IsAppended();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsLCAAncestor() const {
  auto& node = static_cast<const CRTP&>(*this);
  auto& ancestors = node.GetDAG().GetLCAAncestors();
  return ancestors.find(node) != ancestors.end();
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
  if (node.IsMoveTarget()) {
    // then fitch sets can't have changed, but the fitch sets recorded in
    // tree_'s changed fitch set map relative to this node are meant for this
    // node's new parent!
    return {node.GetMATNode().mutations, std::nullopt};
  } else if (node.IsMoveNew()) {
    // Then fitch set changes are relative to the target node
    auto result = dag.GetChangedFitchSetMap().find(
        std::addressof(dag.GetMoveTarget().GetMATNode()));
    return {dag.GetMoveTarget().GetMATNode().mutations,
            result == dag.GetChangedFitchSetMap().end()
                ? std::nullopt
                : std::make_optional(result->second.Copy())};
  } else {
    auto result = dag.GetChangedFitchSetMap().find(std::addressof(node.GetMATNode()));
    return {node.GetMATNode().mutations,
            result == dag.GetChangedFitchSetMap().end()
                ? std::nullopt
                : std::make_optional(result->second.Copy())};
  }
}

template <typename CRTP, typename Tag>
FitchSet FeatureConstView<HypotheticalNode, CRTP, Tag>::GetFitchSet(
    MutationPosition site) const {
  auto node = static_cast<const CRTP&>(*this).Const();
  auto dag = node.GetDAG();
  Assert(site.value < dag.GetReferenceSequence().size());
  auto [old_fitch_sets, changes] = GetFitchSetParts();
  if (old_fitch_sets.find(static_cast<int>(site.value)) ==
      old_fitch_sets.mutations.end()) {
    // if no fitch set is recorded on the corresponding MAT node, we can use
    // a singleton set containing the base at this site in the parent compact genome
    if (changes) {
      return FitchSet({node.GetSingleParent().GetParent().GetCompactGenome().GetBase(
                           site, dag.GetReferenceSequence()) &
                       (~changes.value().at(site).get_decremented() |
                        changes.value().at(site).get_incremented())});
    } else {
      return FitchSet({node.GetSingleParent().GetParent().GetCompactGenome().GetBase(
          site, dag.GetReferenceSequence())});
    }
  } else if (changes.has_value()) {
    return FitchSet(
        (old_fitch_sets.find(static_cast<int>(site.value))->get_all_major_allele() &
         (~changes.value().at(site).get_decremented())) |
        changes.value().at(site).get_incremented());
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
  } else {
    return node.GetSingleParent().GetParent().GetChangedBaseSites().Copy();
  }
}

template <typename CRTP, typename Tag>
CompactGenome FeatureConstView<HypotheticalNode, CRTP, Tag>::ComputeNewCompactGenome()
    const {
  auto node = static_cast<const CRTP&>(*this).Const();
  Assert(node.GetMATNodeId() != NoId);
  ContiguousSet<MutationPosition> changed_base_sites =
      node.GetSitesWithChangedFitchSets();
  const CompactGenome& old_cg = node.GetOld().GetCompactGenome();
  ContiguousMap<MutationPosition, char> cg_changes;
  if (node.IsMATRoot()) {
    // If this node is the root node of the tree, we don't have to worry
    // about parent bases, we just
    // choose a base from each changed fitch set, preferring the base that was
    // already in that site before the SPR move
    ContiguousSet<MutationPosition> focus_sites = GetSitesWithChangedFitchSets();
    for (auto site : focus_sites) {
      FitchSet site_fitch_set = GetFitchSet(site);
      char oldbase = old_cg.GetBase(site, node.GetDAG().GetReferenceSequence());
      if (not site_fitch_set.find(oldbase)) {
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
      FitchSet site_fitch_set = node.GetFitchSet(site);
      char oldbase = old_cg.GetBase(site, node.GetDAG().GetReferenceSequence());
      char parent_base = node.GetSingleParent().GetParent().GetCompactGenome().GetBase(
          site, node.GetDAG().GetReferenceSequence());
      if (site_fitch_set.find(parent_base)) {
        if (oldbase != parent_base) {
          cg_changes.insert({site, parent_base});
          changed_base_sites.insert(site);
        }
      } else if (not site_fitch_set.find(oldbase)) {
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
  if (not IsLCAAncestor()) {
    return (node.GetOld().GetCompactGenome() == node.GetCompactGenome());
    // TODO and node.GetOld().GetLeafSet() == node.GetLeafSet());
  } else {
    return false;
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<HypotheticalNode, CRTP, Tag>::PreorderComputeCompactGenome(
    std::vector<NodeId>& result) const {
  auto& node = static_cast<const CRTP&>(*this);
  if (not node.IsRoot() and not node.IsMoveNew()) {  // TODO
    node.template SetOverlay<Deduplicate<CompactGenome>>();
    node = node.ComputeNewCompactGenome();
  }
  result.push_back(node);
  // If we've reached an anchor node, there's no need to continue down this
  // branch.
  if (not node.IsNonrootAnchorNode()) {
    for (auto child : node.GetChildren()) {
      child.GetChild().PreorderComputeCompactGenome(result);
    }
  }
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetMoveLCA() const {
  auto& self = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetNodeFromMAT(self.data_->move_.LCA->node_id);
};

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetMoveSource() const {
  auto& self = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetNodeFromMAT(self.data_->move_.src->node_id);
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetMoveTarget() const {
  auto& self = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetNodeFromMAT(self.data_->move_.dst->node_id);
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetOldSourceParent() const {
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetMoveSource().GetOld().GetSingleParent().GetParent();
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetOldestChangedNode() const {
  auto& dag = static_cast<const CRTP&>(*this);
  const MAT::Node* oldest_changed_node = std::addressof(dag.GetMoveLCA().GetMATNode());
  for (auto& node_change : dag.GetChangedFitchSetMap()) {
    if (node_change.first->dfs_index < oldest_changed_node->dfs_index) {
      oldest_changed_node = node_change.first;
    }
  }
  return dag.GetNodeFromMAT(oldest_changed_node->node_id);
}

template <typename DAG, typename CRTP, typename Tag>
std::vector<NodeId> FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetFragment()
    const {
  auto& dag = static_cast<const CRTP&>(*this);
  std::vector<NodeId> result;
  auto oldest_changed = dag.GetOldestChangedNode();
  if (oldest_changed.IsRoot()) {
    // we need to add the UA node as the root anchor node of the fragment,
    // somehow
  } else {
    result.push_back(oldest_changed.GetSingleParent().GetParent());
  }
  oldest_changed.PreorderComputeCompactGenome(result);
  return result;
}
template <typename DAG, typename CRTP, typename Tag>
const ContiguousMap<const MAT::Node*,
                    ContiguousMap<MutationPosition, Mutation_Count_Change>>&
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

template <typename DAG, typename CRTP, typename Tag>
void FeatureMutableView<HypotheticalTree<DAG>, CRTP, Tag>::ApplyMove(NodeId src,
                                                                     NodeId dst) const {
  auto& dag = static_cast<const CRTP&>(*this);
  Assert(dag.IsTree());
  auto src_node = dag.Get(src);
  auto dst_node = dag.Get(dst);
  Assert(src_node.GetId() != dst_node.GetId());
  src_node.template SetOverlay<Neighbors>();
  dst_node.template SetOverlay<Neighbors>();
  auto src_parent_edge = src_node.GetSingleParent();
  auto dst_parent_edge = dst_node.GetSingleParent();
  src_parent_edge.template SetOverlay<Endpoints>();
  dst_parent_edge.template SetOverlay<Endpoints>();
  auto src_parent = src_parent_edge.GetParent();
  auto dst_parent = dst_parent_edge.GetParent();
  const bool is_sibling_move = src_parent.GetId() == dst_parent.GetId();
  src_parent.template SetOverlay<Neighbors>();
  if (not is_sibling_move) {
    dst_parent.template SetOverlay<Neighbors>();
  }

  src_parent.RemoveChild(src_parent_edge.GetClade(), src_parent_edge);
  dst_parent.RemoveChild(dst_parent_edge.GetClade(), dst_parent_edge);

  auto new_node = dag.AppendNode();
  auto new_edge = dag.AppendEdge(dst_parent, new_node, dst_parent_edge.GetClade());

  src_parent_edge.Set(new_node, src_node, {0});
  dst_parent_edge.Set(new_node, dst_node, {1});

  new_node.AddEdge(src_parent_edge.GetClade(), src_parent_edge, true);
  new_node.AddEdge(dst_parent_edge.GetClade(), dst_parent_edge, true);
  new_node.AddEdge(new_edge.GetClade(), new_edge, false);
}

template <typename DAG, typename CRTP, typename Tag>
void FeatureMutableView<HypotheticalTree<DAG>, CRTP, Tag>::InitHypotheticalTree(
    const Profitable_Moves& move, const std::vector<Node_With_Major_Allele_Set_Change>&
                                      nodes_with_major_allele_set_change) {
  auto& self = GetFeatureStorage(this);
  Assert(not self.data_);
  auto& dag = static_cast<const CRTP&>(*this);
  dag.ApplyMove(dag.GetNodeFromMAT(move.src->node_id),
                dag.GetNodeFromMAT(move.dst->node_id));
  self.data_ = std::make_unique<typename HypotheticalTree<DAG>::Data>(
      typename HypotheticalTree<DAG>::Data{move, nodes_with_major_allele_set_change});
  if (dag.GetMoveLCA().IsRoot()) {
    return;
  }
  NodeId lca = dag.GetMoveLCA();
  do {
    self.data_->lca_ancestors_.insert(lca);
    lca = dag.Get(lca).GetSingleParent().GetParent();
  } while (not dag.Get(lca).IsRoot());
}

template <typename DAG>
HypotheticalTree<DAG>::Data::Data(const Profitable_Moves& move,
                                  const std::vector<Node_With_Major_Allele_Set_Change>&
                                      nodes_with_major_allele_set_change)
    : move_{move} {
  for (auto& node_with_allele_set_change : nodes_with_major_allele_set_change) {
    Assert(node_with_allele_set_change.node != nullptr);
    ContiguousMap<MutationPosition, Mutation_Count_Change> node_map;
    for (auto& mutation_count_change :
         node_with_allele_set_change.major_allele_set_change) {
      MutationPosition pos = {
          static_cast<size_t>(mutation_count_change.get_position())};
      if (pos.value >= 2147483647) {
        continue;
      }
      node_map.insert({pos, mutation_count_change});
    }
    changed_fitch_set_map_.insert(
        {node_with_allele_set_change.node, std::move(node_map)});
  }
}

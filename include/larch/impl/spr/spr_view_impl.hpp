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
auto FeatureConstView<HypotheticalNode, CRTP, Tag>::GetOld() const {
  auto& node = static_cast<const CRTP&>(*this);
  return node.GetDAG().GetOriginal().Get(node.GetId());
}

template <typename CRTP, typename Tag>
const std::set<MutationPosition>&
FeatureConstView<HypotheticalNode, CRTP, Tag>::GetChangedBaseSites() const {
  return GetFeatureStorage(this).changed_base_sites_;
}

template <typename CRTP, typename Tag>
std::set<MutationPosition>
FeatureConstView<HypotheticalNode, CRTP, Tag>::GetSitesWithChangedFitchSets() const {
  auto& node = static_cast<const CRTP&>(*this);
  std::set<MutationPosition> result;
  std::optional<std::map<MutationPosition, Mutation_Count_Change>> fitch_set_map =
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
          std::optional<std::map<MutationPosition, Mutation_Count_Change>>>
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
                : std::make_optional(result->second)};
  } else {
    auto result = dag.GetChangedFitchSetMap().find(std::addressof(node.GetMATNode()));
    return {node.GetMATNode().mutations, result == dag.GetChangedFitchSetMap().end()
                                             ? std::nullopt
                                             : std::make_optional(result->second)};
  }
}

template <typename CRTP, typename Tag>
FitchSet FeatureConstView<HypotheticalNode, CRTP, Tag>::GetFitchSet(
    MutationPosition site) const {
  auto& node = static_cast<const CRTP&>(*this);
  auto dag = node.GetDAG();
  auto [old_fitch_sets, changes] = GetFitchSetParts();
  if (old_fitch_sets.find(static_cast<int>(site.value)) ==
      old_fitch_sets.mutations.end()) {
    // if no fitch set is recorded on the corresponding MAT node, we can use
    // a singleton set containing the base at this site in the parent compact genome
    // TODO: Does this work if IsSourceNode()?
    return FitchSet({node.GetSingleParent().GetParent().GetCompactGenome().GetBase(
        site, dag.GetReferenceSequence())});
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
std::set<MutationPosition>
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
    return node.GetSingleParent().GetParent().GetChangedBaseSites();
  }
}

template <typename CRTP, typename Tag>
CompactGenome FeatureConstView<HypotheticalNode, CRTP, Tag>::ComputeNewCompactGenome()
    const {
  auto& node = static_cast<const CRTP&>(*this);
  std::set<MutationPosition> changed_base_sites = node.GetSitesWithChangedFitchSets();
  const CompactGenome& old_cg = node.GetOld().GetCompactGenome();
  std::map<MutationPosition, char> cg_changes;
  if (node.IsMATRoot()) {
    // If this node is the root node of the tree, we don't have to worry
    // about parent bases, we just
    // choose a base from each changed fitch set, preferring the base that was
    // already in that site before the SPR move
    std::set<MutationPosition> focus_sites = GetSitesWithChangedFitchSets();
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
    std::set<MutationPosition> focus_sites = node.GetSitesWithChangedFitchSets();
    focus_sites.merge(node.GetParentChangedBaseSites());
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

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetMoveSource() const {
  auto& self = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.Get(NodeId{self.data_->move_.src->node_id});
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetMoveTarget() const {
  auto& self = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.Get(NodeId{self.data_->move_.dst->node_id});
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetOldSourceParent() const {
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetMoveSource().GetOld().GetSingleParent().GetParent();
}

template <typename DAG, typename CRTP, typename Tag>
const std::map<const MAT::Node*, std::map<MutationPosition, Mutation_Count_Change>>&
FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetChangedFitchSetMap() const {
  auto& self = GetFeatureStorage(this);
  return self.data_->changed_fitch_set_map_;
}

template <typename DAG, typename CRTP, typename Tag>
void FeatureMutableView<HypotheticalTree<DAG>, CRTP, Tag>::ApplyMove(NodeId src,
                                                                     NodeId dst) const {
  auto& dag = static_cast<const CRTP&>(*this);
  auto src_node = dag.Get(src);
  auto dst_node = dag.Get(dst);
  auto src_parent_edge = src_node.GetSingleParent();
  auto dst_parent_edge = dst_node.GetSingleParent();
  auto src_parent = src_parent_edge.GetParent();
  auto dst_parent = dst_parent_edge.GetParent();
  auto dst_grandparent_edge = dst_parent.GetSingleParent();

  Assert(src_node.GetId() != dst_node.GetId());
  Assert(src_parent.GetId() != dst_parent.GetId());

  auto new_node = dag.AppendNode();
  auto new_edge = dag.AppendEdge(dst_parent, new_node, {0});

  src_node.Overlay();
  src_parent_edge.Overlay();
  src_parent.Overlay();

  dst_node.Overlay();
  dst_parent_edge.Overlay();
  dst_parent.Overlay();

  src_parent_edge.Set(new_node, src_node, {0});
  dst_parent_edge.Set(new_node, dst_node, {0});

  auto build_connections = [&dag](EdgeId id) {
    auto edge = dag.Get(id);
    edge.GetParent().AddEdge(edge.GetClade(), edge, true);
    edge.GetChild().AddEdge(edge.GetClade(), edge, false);
  };

  for (auto child_edge : src_parent.GetOld().GetChildren()) {
    build_connections(child_edge);
  }

  for (auto child_edge : dst_parent.GetOld().GetChildren()) {
    build_connections(child_edge);
  }

  build_connections(new_edge);
  build_connections(dst_grandparent_edge);
}

template <typename DAG, typename CRTP, typename Tag>
void FeatureMutableView<HypotheticalTree<DAG>, CRTP, Tag>::InitHypotheticalTree(
    const Profitable_Moves& move, const std::vector<Node_With_Major_Allele_Set_Change>&
                                      nodes_with_major_allele_set_change) {
  auto& self = GetFeatureStorage(this);
  self.data_ = std::make_unique<typename HypotheticalTree<DAG>::Data>(
      typename HypotheticalTree<DAG>::Data{move, nodes_with_major_allele_set_change});
  auto& dag = static_cast<const CRTP&>(*this);
  dag.ApplyMove(dag.GetNodeFromMAT(move.src->node_id),
                dag.GetNodeFromMAT(move.dst->node_id));
/*  std::set<size_t> new_nodes;
  for (auto& mat_node : nodes_with_major_allele_set_change) {
    new_nodes.insert(mat_node.node->node_id);
  }
  for (MAT::Node* node = move.src; node != nullptr; node = node->parent) {
    new_nodes.insert(node->node_id);
    if (node == move.LCA) {
      break;
    }
  }
  for (MAT::Node* node = move.dst; node != nullptr; node = node->parent) {
    new_nodes.insert(node->node_id);
    if (node == move.LCA) {
      break;
    }
  }
  for (size_t id : new_nodes) {
    dag.GetMutableNodeFromMAT(id).Overlay();
  }
*/}

template <typename DAG>
HypotheticalTree<DAG>::Data::Data(const Profitable_Moves& move,
                                  const std::vector<Node_With_Major_Allele_Set_Change>&
                                      nodes_with_major_allele_set_change)
    : move_{move} {
  for (auto& node_with_allele_set_change : nodes_with_major_allele_set_change) {
    std::map<MutationPosition, Mutation_Count_Change> node_map;
    for (auto& mutation_count_change :
         node_with_allele_set_change.major_allele_set_change) {
      node_map.insert({{static_cast<size_t>(mutation_count_change.get_position())},
                       mutation_count_change});
    }
    changed_fitch_set_map_.insert(
        {node_with_allele_set_change.node, std::move(node_map)});
  }
}

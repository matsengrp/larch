
template <typename CRTP, typename Tag>
const MAT::Node& FeatureConstView<HypotheticalNode, CRTP, Tag>::GetMATNode() const {
  auto& node = static_cast<const CRTP&>(*this);
  return *node.GetDAG().GetMAT().get_node(node.GetId().value);
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsSource() const {
  auto& node = static_cast<const CRTP&>(*this);
  return node.GetDAG().GetSource().GetId() == node.GetId();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsTarget() const {
  auto& node = static_cast<const CRTP&>(*this);
  return node.GetDAG().GetTarget().GetId() == node.GetId();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<HypotheticalNode, CRTP, Tag>::IsNew() const {
  auto& node = static_cast<const CRTP&>(*this);
  return node.IsOverlaid();
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
  if (node.IsTarget()) {
    // then fitch sets can't have changed, but the fitch sets recorded in
    // tree_'s changed fitch set map relative to this node are meant for this
    // node's new parent!
    return {node.GetMATNode().mutations, std::nullopt};
  } else if (node.IsNew()) {
    // Then fitch set changes are relative to the target node
    auto result =
        dag.GetChangedFitchSetMap().find(std::addressof(dag.GetTarget().GetMATNode()));
    return {dag.GetTarget().GetMATNode().mutations,
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
  auto [old_fitch_sets, changes] = GetFitchSetParts();
  if (old_fitch_sets.find(static_cast<int>(site.value)) ==
      old_fitch_sets.mutations.end()) {
    // if no fitch set is recorded on the corresponding MAT node, we can use
    // a singleton set containing the base at this site in the parent compact genome
    // TODO: Does this work if IsSourceNode()?
    return FitchSet(
        {node.GetSingleParent().GetParent().GetCompactGenome().GetBase(site)});
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
  if (node.IsSourceNode()) {
    // this node used to be below the old source parent, so we need
    // differences between the new parent's new cg, and the old cg of the old
    // parent of the source node.
    const CompactGenome& old_parent_cg = dag.GetOldSourceParent().GetOldCompactGenome();
    // Parent must be the new node:
    const CompactGenome& new_parent_cg = node.GetParent().GetCompactGenome();
    // Imaginary method DifferingSites returns sites at which new_parent_cg
    // and old_parent_cg don't have the same base.
    return old_parent_cg.DifferingSites(new_parent_cg);
  } else {
    return node.GetParent().GetChangedBaseSites();
  }
}

template <typename DAG, typename CRTP, typename Tag>
const MAT::Tree& FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetMAT() const {
  auto& self = GetFeatureStorage(this);
  return self.data_->sample_mat_;
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetSource() const {
  auto& self = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.Get(NodeId{self.data_->move_.src->node_id});
}

template <typename DAG, typename CRTP, typename Tag>
auto FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetTarget() const {
  auto& self = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.Get(NodeId{self.data_->move_.dst->node_id});
}

template <typename DAG, typename CRTP, typename Tag>
const std::map<const MAT::Node*, std::map<MutationPosition, Mutation_Count_Change>>&
FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag>::GetChangedFitchSetMap() const {
  auto& self = GetFeatureStorage(this);
  return self.data_->changed_fitch_set_map_;
}

template <typename DAG, typename CRTP, typename Tag>
void FeatureMutableView<HypotheticalTree<DAG>, CRTP, Tag>::InitHypotheticalTree(
    DAG sample_dag, const MAT::Tree& sample_mat, const Profitable_Moves& move,
    const std::vector<Node_With_Major_Allele_Set_Change>&
        nodes_with_major_allele_set_change) {
  auto& self = GetFeatureStorage(this);
  self.data_ = std::make_unique<typename HypotheticalTree<DAG>::Data>(
      typename HypotheticalTree<DAG>::Data{sample_dag, sample_mat, move,
                                           nodes_with_major_allele_set_change});
}

template <typename DAG>
HypotheticalTree<DAG>::Data::Data(DAG sample_dag, const MAT::Tree& sample_mat,
                                  const Profitable_Moves& move,
                                  const std::vector<Node_With_Major_Allele_Set_Change>&
                                      nodes_with_major_allele_set_change)
    : sample_dag_{sample_dag}, sample_mat_{sample_mat}, move_{move} {
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
#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename CRTP, typename Tag>
NodeId FeatureConstView<MappedNodes, CRTP, Tag>::GetOriginalId() const {
  return GetFeatureStorage(this).original_id_;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<MappedNodes, CRTP, Tag>::SetOriginalId(NodeId id) const {
  GetFeatureStorage(this).original_id_ = id;
  auto& node = static_cast<const CRTP&>(*this);
  node.GetDAG()
      .template GetFeatureExtraStorage<NodeId, MappedNodes>()
      .reverse_map_.Insert(id, node.GetId());
}

template <typename CRTP>
auto ExtraFeatureConstView<MappedNodes, CRTP>::GetMappedNode(NodeId original_id) const {
  auto& dag = static_cast<const CRTP&>(*this);
  NodeId id =
      dag.template GetFeatureExtraStorage<NodeId, MappedNodes>().reverse_map_.at(
          original_id);
  return dag.Get(id);
}

template <typename CRTP>
auto ExtraFeatureMutableView<MappedNodes, CRTP>::GetMutableMappedNode(
    NodeId original_id) const {
  auto& dag = static_cast<const CRTP&>(*this);
  NodeId id =
      dag.template GetFeatureExtraStorage<NodeId, MappedNodes>().reverse_map_.at(
          original_id);
  return dag.Get(id);
}
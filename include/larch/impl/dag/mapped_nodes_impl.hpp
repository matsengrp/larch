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
}
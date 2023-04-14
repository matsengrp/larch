#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename CRTP, typename Tag>
auto FeatureConstView<Endpoints, CRTP, Tag>::GetParent() const {
  auto dag = static_cast<const CRTP&>(*this).GetDAG();
  return typename decltype(dag)::NodeView{dag, GetParentId()};
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Endpoints, CRTP, Tag>::GetChild() const {
  auto dag = static_cast<const CRTP&>(*this).GetDAG();
  return typename decltype(dag)::NodeView{dag, GetChildId()};
}

template <typename CRTP, typename Tag>
CladeIdx FeatureConstView<Endpoints, CRTP, Tag>::GetClade() const {
  return GetFeatureStorage(this).clade_;
}

template <typename CRTP, typename Tag>
NodeId FeatureConstView<Endpoints, CRTP, Tag>::GetParentId() const {
  return GetFeatureStorage(this).parent_;
}

template <typename CRTP, typename Tag>
NodeId FeatureConstView<Endpoints, CRTP, Tag>::GetChildId() const {
  return GetFeatureStorage(this).child_;
}

template <typename CRTP, typename Tag>
std::pair<NodeId, NodeId> FeatureConstView<Endpoints, CRTP, Tag>::GetNodeIds() const {
  return {GetParentId(), GetChildId()};
}

template <typename CRTP, typename Tag>
bool FeatureConstView<Endpoints, CRTP, Tag>::IsUA() const {
  return GetParent().IsUA();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<Endpoints, CRTP, Tag>::IsLeaf() const {
  return GetChild().IsLeaf();
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Endpoints, CRTP, Tag>::Set(NodeId parent, NodeId child,
                                                   CladeIdx clade) const {
  Assert(parent.value != NoId);
  Assert(child.value != NoId);
  Assert(parent.value != child.value);
  Assert(clade.value != NoId);
  GetFeatureStorage(this).parent_ = parent;
  GetFeatureStorage(this).child_ = child;
  GetFeatureStorage(this).clade_ = clade;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Endpoints, CRTP, Tag>::SetParent(NodeId parent) const {
  Assert(parent.value != NoId);
  GetFeatureStorage(this).parent_ = parent;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Endpoints, CRTP, Tag>::SetChild(NodeId child) const {
  Assert(child.value != NoId);
  GetFeatureStorage(this).child_ = child;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Endpoints, CRTP, Tag>::SetClade(CladeIdx clade) const {
  Assert(clade.value != NoId);
  GetFeatureStorage(this).clade_ = clade;
}

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
  return GetStorageClade();
}

template <typename CRTP, typename Tag>
NodeId FeatureConstView<Endpoints, CRTP, Tag>::GetParentId() const {
  return GetStorageParent();
}

template <typename CRTP, typename Tag>
NodeId FeatureConstView<Endpoints, CRTP, Tag>::GetChildId() const {
  return GetStorageChild();
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
bool FeatureConstView<Endpoints, CRTP, Tag>::IsTreeRoot() const {
  return GetParent().IsTreeRoot();
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
  GetStorage().SetParent(static_cast<const CRTP*>(this), parent);
  GetStorage().SetChild(static_cast<const CRTP*>(this), child);
  GetStorage().SetClade(static_cast<const CRTP*>(this), clade);
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Endpoints, CRTP, Tag>::SetParent(NodeId parent) const {
  Assert(parent.value != NoId);
  GetStorage().SetParent(static_cast<const CRTP*>(this), parent);
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Endpoints, CRTP, Tag>::SetChild(NodeId child) const {
  Assert(child.value != NoId);
  GetStorage().SetChild(static_cast<const CRTP*>(this), child);
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Endpoints, CRTP, Tag>::SetClade(CladeIdx clade) const {
  Assert(clade.value != NoId);
  GetStorage().SetClade(static_cast<const CRTP*>(this), clade);
}

#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename DAG>
EdgeView<DAG>::EdgeView(DAG dag, EdgeId id) : dag_{dag}, id_{id} {}

template <typename DAG>
EdgeView<DAG>::operator EdgeId() {
  return id_;
}

template <typename DAG>
EdgeView<DAG>::operator CladeIdx() {
  return GetStorage().clade_;
}

template <typename DAG>
auto& EdgeView<DAG>::GetDAG() {
  return dag_;
}

template <typename DAG>
EdgeId EdgeView<DAG>::GetId() {
  return id_;
}

template <typename DAG>
auto EdgeView<DAG>::GetParent() {
  return Node{dag_, GetParentId()};
}

template <typename DAG>
auto EdgeView<DAG>::GetChild() {
  return Node{dag_, GetChildId()};
}

template <typename DAG>
CladeIdx EdgeView<DAG>::GetClade() {
  return GetStorage().clade_;
}

template <typename DAG>
NodeId EdgeView<DAG>::GetParentId() {
  return GetStorage().parent_;
}

template <typename DAG>
NodeId EdgeView<DAG>::GetChildId() {
  return GetStorage().child_;
}

template <typename DAG>
std::pair<NodeId, NodeId> EdgeView<DAG>::GetNodeIds() {
  return {GetStorage().parent_, GetStorage().child_};
}

template <typename DAG>
bool EdgeView<DAG>::IsRoot() {
  return GetParent().IsRoot();
}

template <typename DAG>
bool EdgeView<DAG>::IsLeaf() {
  return GetChild().IsLeaf();
}

template <typename DAG>
auto& EdgeView<DAG>::GetStorage() const {
  return dag_.storage_.edges_.edges_.at(id_.value);
}

template <typename DAG>
template <typename Feature>
auto& EdgeView<DAG>::GetFeatureStorage() const {
  if constexpr (DAG::StorageType::EdgesContainerType::template contains_feature<Feature>) {
    return dag_.storage_.edges_.template GetFeatureAt<Feature>(id_);
  } else {
    return std::get<Feature>(dag_.storage_.edges_.edges_.at(id_.value).features_);
  }
}

template <std::size_t Index, typename DAG>
std::tuple_element_t<Index, EdgeView<DAG>> get(EdgeView<DAG> edge) {
  if constexpr (Index == 0) {
    return edge.GetParent();
  }
  if constexpr (Index == 1) {
    return edge.GetChild();
  }
}

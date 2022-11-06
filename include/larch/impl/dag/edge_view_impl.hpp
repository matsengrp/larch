#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename DAGType, typename... Features>
EdgeView<DAGType, Features...>::EdgeView(DAGType dag, EdgeId id) : dag_{dag}, id_{id} {}

template <typename DAGType, typename... Features>
EdgeView<DAGType, Features...>::operator EdgeId() const {
  return id_;
}

template <typename DAGType, typename... Features>
EdgeView<DAGType, Features...>::operator CladeIdx() const {
  return GetStorage().clade_;
}

template <typename DAGType, typename... Features>
auto& EdgeView<DAGType, Features...>::GetDAG() const {
  return dag_;
}

template <typename DAGType, typename... Features>
EdgeId EdgeView<DAGType, Features...>::GetId() const {
  return id_;
}

template <typename DAGType, typename... Features>
auto EdgeView<DAGType, Features...>::GetParent() const {
  return Node{dag_, GetParentId()};
}

template <typename DAGType, typename... Features>
auto EdgeView<DAGType, Features...>::GetChild() const {
  return Node{dag_, GetChildId()};
}

template <typename DAGType, typename... Features>
CladeIdx EdgeView<DAGType, Features...>::GetClade() const {
  return GetStorage().clade_;
}

template <typename DAGType, typename... Features>
NodeId EdgeView<DAGType, Features...>::GetParentId() const {
  return GetStorage().parent_;
}

template <typename DAGType, typename... Features>
NodeId EdgeView<DAGType, Features...>::GetChildId() const {
  return GetStorage().child_;
}

template <typename DAGType, typename... Features>
std::pair<NodeId, NodeId> EdgeView<DAGType, Features...>::GetNodeIds() const {
  return {GetStorage().parent_, GetStorage().child_};
}

// template <typename DAGType, typename... Features>
// bool EdgeView<DAGType, Features...>::IsRoot() const;

// template <typename DAGType, typename... Features>
// bool EdgeView<DAGType, Features...>::IsLeaf() const;

template <typename DAGType, typename... Features>
auto& EdgeView<DAGType, Features...>::GetStorage() const {
  return dag_.storage_.edges_.edges_.at(id_.value);
}

template <std::size_t Index, typename DAGType, typename... Features>
std::tuple_element_t<Index, EdgeView<DAGType, Features...>> get(
    EdgeView<DAGType, Features...> edge) {
  if constexpr (Index == 0) {
    return edge.GetParent();
  }
  if constexpr (Index == 1) {
    return edge.GetChild();
  }
}

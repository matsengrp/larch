template <typename T>
TraverseValue<T>::TraverseValue(T dag, NodeId node, EdgeId edge)
    : dag_{dag}, node_{node}, edge_{edge} {}

template <typename T>
NodeView<T> TraverseValue<T>::GetNode() const {
  return dag_.Get(node_);
}

template <typename T>
EdgeView<T> TraverseValue<T>::GetEdge() const {
  return dag_.Get(edge_);
}

template <typename T>
TraverseValue<T>::operator MutableNode() const {
  return GetNode();
}

template <typename T>
TraverseValue<T>::operator MutableEdge() const {
  return GetEdge();
}

template <typename T>
TraverseValue<T>::operator Node() const {
  return GetNode();
}

template <typename T>
TraverseValue<T>::operator Edge() const {
  return GetEdge();
}

template <typename T>
bool TraverseValue<T>::IsRoot() const {
  return GetNode().IsRoot();
}

template <typename T>
bool TraverseValue<T>::IsLeaf() const {
  return GetNode().IsLeaf();
}

namespace std {
template <typename T>
struct tuple_size<::TraverseValue<T>> : integral_constant<size_t, 2> {};

template <size_t Index, typename T>
struct tuple_element<Index, ::TraverseValue<T>>
    : tuple_element<Index, tuple<NodeView<T>, EdgeView<T>>> {};
}  // namespace std

template <std::size_t Index, typename T>
std::tuple_element_t<Index, TraverseValue<T>> get(TraverseValue<T> tv) {
  if constexpr (Index == 0) return tv.GetNode();
  if constexpr (Index == 1) return tv.GetEdge();
}

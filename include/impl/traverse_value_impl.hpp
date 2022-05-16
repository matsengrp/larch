template <typename T>
TraverseValue<T>::TraverseValue(T dag, NodeId node, EdgeId edge) :
		dag_{dag}, node_{node}, edge_{edge} {}

template <typename T>
NodeView<T> TraverseValue<T>::GetNode() const { return dag_.GetNode(node_); }

template <typename T>
EdgeView<T> TraverseValue<T>::GetEdge() const { return dag_.GetEdge(edge_); }

template <typename T>
TraverseValue<T>::operator MutableNode() const { return GetNode(); }

template <typename T>
TraverseValue<T>::operator MutableEdge() const { return GetEdge(); }

template <typename T>
TraverseValue<T>::operator Node() const { return GetNode(); }

template <typename T>
TraverseValue<T>::operator Edge() const { return GetEdge(); }

template <typename T>
bool TraverseValue<T>::IsRoot() const { return GetNode().IsRoot(); }

template <typename T>
bool TraverseValue<T>::IsLeaf() const { return GetNode().IsLeaf(); }

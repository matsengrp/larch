#include "edge.hpp"

#include "dag.hpp"

template <typename T>
EdgeView<T>::EdgeView(T dag, EdgeId id) : dag_{dag}, id_{id} {
  static_assert(std::is_same_v<T, DAG&> or std::is_same_v<T, const DAG&>);
  Assert(id.value != NoId);
  Assert(id.value < dag_.edges_.size());
}

template <typename T>
EdgeView<T>::operator Edge() const {
  return {dag_, id_};
}

template <typename T>
EdgeView<T>::operator EdgeId() const {
  return id_;
}

template <typename T>
EdgeView<T>::operator CladeIdx() const {
  return GetClade();
}

template <typename T>
T EdgeView<T>::GetDAG() const {
  return dag_;
}

template <typename T>
EdgeId EdgeView<T>::GetId() const {
  return id_;
}

template <typename T>
typename EdgeView<T>::NodeType EdgeView<T>::GetParent() const {
  return {dag_, GetStorage().GetParent()};
}

template <typename T>
typename EdgeView<T>::NodeType EdgeView<T>::GetChild() const {
  return {dag_, GetStorage().GetChild()};
}

template <typename T>
CladeIdx EdgeView<T>::GetClade() const {
  return GetStorage().GetClade();
}

template <typename T>
NodeId EdgeView<T>::GetParentId() const {
  return GetParent();
}

template <typename T>
NodeId EdgeView<T>::GetChildId() const {
  return GetChild();
}

template <typename T>
std::pair<NodeId, NodeId> EdgeView<T>::GetNodeIds() const {
  return {GetParentId(), GetChildId()};
}

template <typename T>
bool EdgeView<T>::IsRoot() const {
  return GetParent().IsRoot();
}

template <typename T>
bool EdgeView<T>::IsLeaf() const {
  return GetChild().IsLeaf();
}

template <typename T>
std::optional<EdgeView<T>> EdgeView<T>::FindNextSibling() const {
  auto [parent, child] = *this;
  auto children = parent.GetChildren();
  for (auto i = children.begin(); i != children.end(); ++i) {
    if ((*i).GetChildId() == child.GetId()) {
      if (++i != children.end()) {
        return *i;
      }
      break;
    }
  }
  return std::nullopt;
}

template <typename T>
const auto& EdgeView<T>::GetStorage() const {
  return dag_.edges_.at(id_.value);
}

template <typename T>
auto& EdgeView<T>::GetStorage() {
  return dag_.edges_.at(id_.value);
}

template class EdgeView<DAG&>;
template class EdgeView<const DAG&>;

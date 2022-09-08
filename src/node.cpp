#include "node.hpp"

#include "dag.hpp"

template <typename T>
NodeView<T>::NodeView(T dag, NodeId id) : dag_{dag}, id_{id} {
  static_assert(std::is_same_v<T, DAG&> or std::is_same_v<T, const DAG&>);
  Assert(id.value != NoId);
  Assert(id.value < dag_.nodes_.size());
}

template <typename T>
NodeView<T>::operator Node() const {
  return {dag_, id_};
}

template <typename T>
NodeView<T>::operator NodeId() const {
  return id_;
}

template <typename T>
T NodeView<T>::GetDAG() const {
  return dag_;
}

template <typename T>
NodeId NodeView<T>::GetId() const {
  return id_;
}

template <typename T>
typename NodeView<T>::EdgeType NodeView<T>::GetSingleParent() const {
  Assert(GetParents().size() == 1);
  return *GetParents().begin();
}

template <typename T>
bool NodeView<T>::IsRoot() const {
  return GetStorage().GetParents().empty();
}

template <typename T>
bool NodeView<T>::IsLeaf() const {
  auto children = GetChildren();
  return children.begin() == children.end();
}

template <typename T>
void NodeView<T>::AddParentEdge(Edge edge) const {
  if constexpr (is_mutable) {
    GetStorage().AddEdge(edge.GetClade(), edge.GetId(), false);
  }
}

template <typename T>
void NodeView<T>::AddChildEdge(Edge edge) const {
  if constexpr (is_mutable) {
    GetStorage().AddEdge(edge.GetClade(), edge.GetId(), true);
  }
}

template <typename T>
void NodeView<T>::RemoveParentEdge(Edge edge) const {
  if constexpr (is_mutable) {
    GetStorage().RemoveEdge(edge, false);
  }
}

template <typename T>
const std::optional<std::string> NodeView<T>::GetSampleId() const {
    return GetStorage().GetSampleId();
}

template <typename T>
void NodeView<T>::SetSampleId(std::optional<std::string> sample_id) {
  if constexpr (is_mutable) {
    GetStorage().SetSampleId(sample_id);
  }
}

template class NodeView<DAG&>;
template class NodeView<const DAG&>;

#include "history_dag_node_storage.hpp"

template <typename T>
NodeView<T>::NodeView(T dag, NodeId id) : dag_{dag}, id_{id} {
  static_assert(std::is_same_v<T, HistoryDAG&> or std::is_same_v<T, const HistoryDAG&>);
  assert(id.value != NoId);
  assert(id.value < dag_.nodes_.size());
}

template <typename T>
NodeView<T>::operator Node() const {
  return {dag_, id_};
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
auto NodeView<T>::GetParents() const {
  return GetStorage().parents_ | ranges::views::transform([this](EdgeId idx) {
           return EdgeType{dag_, idx};
         });
}

template <typename T>
auto NodeView<T>::GetClades() const {
  return GetStorage().clades_ |
         ranges::views::transform([this](const std::vector<EdgeId>& clade) {
           return clade | ranges::views::transform([this](EdgeId idx) {
                    return EdgeType{dag_, idx};
                  });
         });
}

template <typename T>
auto NodeView<T>::GetChildren() const {
  return GetClades() | ranges::views::join;
}

template <typename T>
typename NodeView<T>::EdgeType NodeView<T>::GetSingleParent() const {
  assert(GetParents().size() == 1);
  return *GetParents().begin();
}

template <typename T>
typename NodeView<T>::EdgeType NodeView<T>::GetFirstParent() const {
  assert(GetParents().size() > 0);
  return *GetParents().begin();
}

template <typename T>
bool NodeView<T>::IsRoot() const {
  return GetParents().empty();
}

template <typename T>
bool NodeView<T>::IsLeaf() const {
  return GetChildren().begin() == GetChildren().end();
}

template <typename T>
void NodeView<T>::AddParentEdge(Edge edge) {
  GetStorage().AddEdge(edge.GetClade(), edge.GetId(), false);
}

template <typename T>
void NodeView<T>::AddChildEdge(Edge edge) {
  GetStorage().AddEdge(edge.GetClade(), edge.GetId(), true);
}

template <typename T>
void NodeView<T>::RemoveParentEdge(Edge edge) {
  GetStorage().RemoveEdge(edge, false);
}

template <typename T>
auto& NodeView<T>::GetStorage() const {
  return dag_.nodes_.at(id_.value);
}

template <typename T>
inline bool operator==(NodeView<T> lhs, NodeView<T> rhs) {
  return std::addressof(lhs.dag_) == std::addressof(rhs.dag_) && lhs.id_ == rhs.id_;
}

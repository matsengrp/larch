#include <range/v3/view/join.hpp>

// Return a range containing this node's parent Edges
template <typename T>
auto NodeView<T>::GetParents() const {
  return GetStorage().GetParents() | Transform::ToEdges(dag_);
}

// Return a range containing this node's child Edges, grouped in sub-ranges by clade
template <typename T>
auto NodeView<T>::GetClades() const {
  return GetStorage().GetClades() |
         ranges::views::transform([this](const std::vector<EdgeId>& clade) {
           return clade | Transform::ToEdges(dag_);
         });
}

// Return a range containing this node's child Edges
template <typename T>
auto NodeView<T>::GetChildren() const {
  return GetClades() | ranges::views::join;
}

template <typename T>
auto& NodeView<T>::GetStorage() const {
  return dag_.nodes_.at(id_.value);
}

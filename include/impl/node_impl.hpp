#include <range/v3/view/join.hpp>

template <typename T>
auto NodeView<T>::GetParents() const {
  return GetStorage().parents_ | Transform::ToEdges(dag_);
}

template <typename T>
auto NodeView<T>::GetClades() const {
  return GetStorage().clades_ |
         ranges::views::transform([this](const std::vector<EdgeId>& clade) {
           return clade | Transform::ToEdges(dag_);
         });
}

template <typename T>
auto NodeView<T>::GetChildren() const {
  return GetClades() | ranges::views::join;
}

template <typename T>
auto& NodeView<T>::GetStorage() const {
  return dag_.nodes_.at(id_.value);
}

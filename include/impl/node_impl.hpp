// Functions defined here are documented where declared in `include/node.hpp`
#include <range/v3/view/join.hpp>

template <typename T>
auto NodeView<T>::GetParents() const {
  return GetStorage().GetParents() | Transform::ToEdges(dag_);
}

template <typename T>
auto NodeView<T>::GetClades() const {
  return GetStorage().GetClades() |
         ranges::views::transform([this](const std::vector<EdgeId>& clade) {
           return clade | Transform::ToEdges(dag_);
         });
}

template <typename T>
auto NodeView<T>::GetClade(CladeIdx clade) const {
  return GetStorage().GetClades().at(clade.value) | Transform::ToEdges(dag_);
}

template <typename T>
size_t NodeView<T>::GetCladesCount() const {
  return GetStorage().GetClades().size();
}

template <typename T>
auto NodeView<T>::GetChildren() const {
  return GetClades() | ranges::views::join;
}

template <typename T>
auto& NodeView<T>::GetStorage() const {
  return dag_.nodes_.at(id_.value);
}

#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Storage, typename... Features>
template <typename Feature>
Feature& DefaultEdgesContainer<Storage, Features...>::GetFeatureAt(EdgeId id) {
  return std::get<std::vector<Feature>>(features_).at(id.value);
}

template <typename Storage, typename... Features>
template <typename Feature>
const Feature& DefaultEdgesContainer<Storage, Features...>::GetFeatureAt(
    EdgeId id) const {
  return std::get<std::vector<Feature>>(features_).at(id.value);
}

template <typename Storage, typename... Features>
Storage& DefaultEdgesContainer<Storage, Features...>::AddEdge(EdgeId id) {
  return GetOrInsert(edges_, id);
}

template <typename Storage, typename... Features>
auto DefaultEdgesContainer<Storage, Features...>::View() {
  return edges_ | ranges::views::all;
}

template <typename Storage, typename... Features>
auto DefaultEdgesContainer<Storage, Features...>::View() const {
  return edges_ | ranges::views::all;
}

template <typename Storage, typename... Features>
size_t DefaultEdgesContainer<Storage, Features...>::Count() const {
  return edges_.size();
}

template <typename Storage, typename... Features>
Storage& DefaultEdgesContainer<Storage, Features...>::EdgeAt(EdgeId id) {
  return edges_.at(id.value);
}

template <typename Storage, typename... Features>
const Storage& DefaultEdgesContainer<Storage, Features...>::EdgeAt(EdgeId id) const {
  return edges_.at(id.value);
}
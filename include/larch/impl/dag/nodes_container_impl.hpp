#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Storage, typename... Features>
template <typename Feature, typename Id>
Feature& DefaultNodesContainer<Storage, Features...>::GetFeatureAt(Id id) {
  return std::get<std::vector<Feature>>(features_).at(id.value);
}

template <typename Storage, typename... Features>
template <typename Feature, typename Id>
const Feature& DefaultNodesContainer<Storage, Features...>::GetFeatureAt(Id id) const {
  return std::get<std::vector<Feature>>(features_).at(id.value);
}
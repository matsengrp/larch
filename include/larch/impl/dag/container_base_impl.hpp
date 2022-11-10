#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Id, template <typename> typename ViewType, typename Storage,
          typename... Features>
template <typename Feature>
Feature& DefaultContainerBase<Id, ViewType, Storage, Features...>::GetFeatureAt(Id id) {
  return std::get<std::vector<Feature>>(features_).at(id.value);
}

template <typename Id, template <typename> typename ViewType, typename Storage,
          typename... Features>
template <typename Feature>
const Feature& DefaultContainerBase<Id, ViewType, Storage, Features...>::GetFeatureAt(
    Id id) const {
  return std::get<std::vector<Feature>>(features_).at(id.value);
}

template <typename Id, template <typename> typename ViewType, typename Storage,
          typename... Features>
template <typename Feature>
typename Feature::GlobalData&
DefaultContainerBase<Id, ViewType, Storage, Features...>::GetGlobalData() {
  return std::get<Feature::GlobalData>(global_data_);
}

template <typename Id, template <typename> typename ViewType, typename Storage,
          typename... Features>
template <typename Feature>
const typename Feature::GlobalData&
DefaultContainerBase<Id, ViewType, Storage, Features...>::GetGlobalData() const {
  return std::get<Feature::GlobalData>(global_data_);
}
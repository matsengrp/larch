#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Id, template <typename> typename ViewType, typename Storage,
          typename... Features>
template <typename Feature>
const Feature& DefaultContainerBase<Id, ViewType, Storage, Features...>::GetFeatureAt(
    Id id) const {
  if constexpr (contains_global_data<Feature>) {
    return std::get<Feature>(features_global_data_)
        .Get(std::get<std::vector<Feature>>(features_).at(id.value), id);
  } else {
    return std::get<std::vector<Feature>>(features_).at(id.value);
  }
}

template <typename Id, template <typename> typename ViewType, typename Storage,
          typename... Features>
template <typename Feature>
void DefaultContainerBase<Id, ViewType, Storage, Features...>::SetFeatureAt(
    Id id, Feature&& feature) {
  if constexpr (contains_global_data<Feature>) {
    std::get<Feature>(features_global_data_)
        .Set(std::get<std::vector<Feature>>(features_).at(id.value), id,
             std::forward<Feature>(feature));
  } else {
    std::get<std::vector<Feature>>(features_).at(id.value) =
        std::forward<Feature>(feature);
  }
}

template <typename Id, template <typename> typename ViewType, typename Storage,
          typename... Features>
template <typename Feature>
auto& DefaultContainerBase<Id, ViewType, Storage, Features...>::GetFeatureGlobalData() {
  using Type =
      typename FeatureTraits<Feature,
                             ViewType<DAGView<Storage>>>::template GlobalData<Id>;
  return std::get<Type>(features_global_data_);
}

template <typename Id, template <typename> typename ViewType, typename Storage,
          typename... Features>
template <typename Feature>
const auto&
DefaultContainerBase<Id, ViewType, Storage, Features...>::GetFeatureGlobalData() const {
  using Type =
      typename FeatureTraits<Feature,
                             ViewType<DAGView<Storage>>>::template GlobalData<Id>;
  return std::get<Type>(features_global_data_);
}
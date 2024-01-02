#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
template <typename Feature>
inline constexpr bool ElementsContainer<C, ElementStorageT, IdCont,
                                        Features...>::contains_element_feature =
    tuple_contains_v<std::tuple<Features...>, Feature> or
    ElementStorageT::template contains_element_feature<Feature>;

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
size_t ElementsContainer<C, ElementStorageT, IdCont, Features...>::GetCount() const {
  return elements_storage_.size();
}

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
Id<C> ElementsContainer<C, ElementStorageT, IdCont, Features...>::Append() {
  Id<C> result{GetCount()};
  elements_storage_.push_back({});
  features_storage_.push_back({});
  return result;
}

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
void ElementsContainer<C, ElementStorageT, IdCont, Features...>::Add(Id<C> id) {
  elements_storage_[id];
  features_storage_[id];
}

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
void ElementsContainer<C, ElementStorageT, IdCont, Features...>::Initialize(
    size_t size) {
  elements_storage_.resize(size);
  features_storage_.resize(size);
}

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
void ElementsContainer<C, ElementStorageT, IdCont, Features...>::Clear() {
  elements_storage_.clear();
  features_storage_.clear();
}

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
template <typename Feature>
auto& ElementsContainer<C, ElementStorageT, IdCont, Features...>::GetFeatureStorage(
    Id<C> id) {
  if constexpr (tuple_contains_v<std::tuple<Features...>, Feature>) {
    return std::get<Feature>(features_storage_.at(id));
  } else {
    return elements_storage_.at(id).template GetFeatureStorage<Feature>();
  }
}

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
template <typename Feature>
const auto& ElementsContainer<C, ElementStorageT, IdCont,
                              Features...>::GetFeatureStorage(Id<C> id) const {
  if constexpr (tuple_contains_v<std::tuple<Features...>, Feature>) {
    return std::get<Feature>(features_storage_.at(id));
  } else {
    return elements_storage_.at(id).template GetFeatureStorage<Feature>();
  }
}

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
template <typename Feature>
auto& ElementsContainer<C, ElementStorageT, IdCont,
                        Features...>::GetFeatureExtraStorage() {
  if constexpr (tuple_contains_v<decltype(extra_features_storage_), Feature>) {
    return std::get<ExtraFeatureStorage<Feature>>(extra_features_storage_);
  } else {
    return std::get<ExtraFeatureStorage<Feature>>(elements_extra_features_storage_);
  }
}

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
template <typename Feature>
const auto& ElementsContainer<C, ElementStorageT, IdCont,
                              Features...>::GetFeatureExtraStorage() const {
  if constexpr (tuple_contains_v<decltype(extra_features_storage_), Feature>) {
    return std::get<ExtraFeatureStorage<Feature>>(extra_features_storage_);
  } else {
    return std::get<ExtraFeatureStorage<Feature>>(elements_extra_features_storage_);
  }
}

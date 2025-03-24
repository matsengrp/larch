#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
template <typename Feature>
inline constexpr bool ElementsContainer<C, ElementStorageT, IdCont,
                                        Features...>::contains_element_feature =
    tuple_contains_v<std::tuple<Features...>, Feature, FeatureEquivalent> or
    ElementStorageT::template contains_element_feature<Feature>;

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
template <typename VT>
size_t ElementsContainer<C, ElementStorageT, IdCont, Features...>::GetCount() const {
  return elements_storage_.size();
}

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
Id<C> ElementsContainer<C, ElementStorageT, IdCont, Features...>::Append() {
  Id<C> result{GetCount<void>()};
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
template <typename Feature, typename E>
auto ElementsContainer<C, ElementStorageT, IdCont, Features...>::GetFeatureStorage(
    Id<C> id, E) {
  if constexpr (tuple_contains_v<std::tuple<Features...>, Feature, FeatureEquivalent>) {
    return std::ref(tuple_get<Feature, FeatureEquivalent>(features_storage_.at(id)));
  } else {
    return elements_storage_.at(id).template GetFeatureStorage<Feature>();
  }
}

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
template <typename Feature, typename E>
auto ElementsContainer<C, ElementStorageT, IdCont, Features...>::GetFeatureStorage(
    Id<C> id, E) const {
  if constexpr (tuple_contains_v<std::tuple<Features...>, Feature, FeatureEquivalent>) {
    return std::cref(tuple_get<Feature, FeatureEquivalent>(features_storage_.at(id)));
  } else {
    return elements_storage_.at(id).template GetFeatureStorage<Feature>();
  }
}

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
template <typename Feature>
auto ElementsContainer<C, ElementStorageT, IdCont,
                       Features...>::GetFeatureExtraStorage() {
  if constexpr (tuple_contains_v<decltype(extra_features_storage_), Feature,
                                 FeatureEquivalent>) {
    return std::ref(tuple_get<ExtraFeatureStorage<Feature>, FeatureEquivalent>(
        extra_features_storage_));
  } else {
    return std::ref(tuple_get<ExtraFeatureStorage<Feature>, FeatureEquivalent>(
        elements_extra_features_storage_));
  }
}

template <Component C, typename ElementStorageT, IdContinuity IdCont,
          typename... Features>
template <typename Feature>
auto ElementsContainer<C, ElementStorageT, IdCont,
                       Features...>::GetFeatureExtraStorage() const {
  if constexpr (tuple_contains_v<decltype(extra_features_storage_), Feature,
                                 FeatureEquivalent>) {
    return std::cref(tuple_get<ExtraFeatureStorage<Feature>, FeatureEquivalent>(
        extra_features_storage_));
  } else {
    return std::cref(tuple_get<ExtraFeatureStorage<Feature>, FeatureEquivalent>(
        elements_extra_features_storage_));
  }
}

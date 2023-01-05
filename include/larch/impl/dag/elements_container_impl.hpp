#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Id, typename ElementStorageT, typename... Features>
template <typename Feature>
inline constexpr bool
    ElementsContainer<Id, ElementStorageT, Features...>::contains_element_feature =
        tuple_contains_v<std::tuple<Features...>, Feature> or
        ElementStorageT::template contains_element_feature<Feature>;

template <typename Id, typename ElementStorageT, typename... Features>
size_t ElementsContainer<Id, ElementStorageT, Features...>::GetCount() const {
  return elements_storage_.size();
}

template <typename Id, typename ElementStorageT, typename... Features>
Id ElementsContainer<Id, ElementStorageT, Features...>::Append() {
  Id result{GetCount()};
  elements_storage_.push_back({});
  features_storage_.push_back({});
  return result;
}

template <typename Id, typename ElementStorageT, typename... Features>
void ElementsContainer<Id, ElementStorageT, Features...>::Add(Id id) {
  std::ignore = GetOrInsert(elements_storage_, id);
  std::ignore = GetOrInsert(features_storage_, id);
}

template <typename Id, typename ElementStorageT, typename... Features>
void ElementsContainer<Id, ElementStorageT, Features...>::Initialize(size_t size) {
  elements_storage_.resize(size);
  features_storage_.resize(size);
}

template <typename Id, typename ElementStorageT, typename... Features>
template <typename Feature>
auto& ElementsContainer<Id, ElementStorageT, Features...>::GetFeatureStorage(Id id) {
  if constexpr (tuple_contains_v<std::tuple<Features...>, Feature>) {
    return std::get<Feature>(features_storage_.at(id.value));
  } else {
    return elements_storage_.at(id.value).template GetFeatureStorage<Feature>();
  }
}

template <typename Id, typename ElementStorageT, typename... Features>
template <typename Feature>
const auto& ElementsContainer<Id, ElementStorageT, Features...>::GetFeatureStorage(
    Id id) const {
  if constexpr (tuple_contains_v<std::tuple<Features...>, Feature>) {
    return std::get<Feature>(features_storage_.at(id.value));
  } else {
    return elements_storage_.at(id.value).template GetFeatureStorage<Feature>();
  }
}

template <typename Id, typename ElementStorageT, typename... Features>
template <typename Feature>
auto& ElementsContainer<Id, ElementStorageT, Features...>::GetFeatureExtraStorage() {
  if constexpr (tuple_contains_v<decltype(extra_features_storage_), Feature>) {
    return std::get<ExtraFeatureStorage<Feature>>(extra_features_storage_);
  } else {
    return std::get<ExtraFeatureStorage<Feature>>(elements_extra_features_storage_);
  }
}

template <typename Id, typename ElementStorageT, typename... Features>
template <typename Feature>
const auto&
ElementsContainer<Id, ElementStorageT, Features...>::GetFeatureExtraStorage() const {
  if constexpr (tuple_contains_v<decltype(extra_features_storage_), Feature>) {
    return std::get<ExtraFeatureStorage<Feature>>(extra_features_storage_);
  } else {
    return std::get<ExtraFeatureStorage<Feature>>(elements_extra_features_storage_);
  }
}

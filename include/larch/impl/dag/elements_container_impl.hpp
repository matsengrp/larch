#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Id, typename ES, typename... Fs>
template <typename Feature>
inline constexpr bool ElementsContainer<Id, ES, Fs...>::contains_element_feature =
    tuple_contains_v<std::tuple<Fs...>, Feature> or
    ES::template contains_element_feature<Feature>;

template <typename Id, typename ES, typename... Fs>
size_t ElementsContainer<Id, ES, Fs...>::GetCount() const {
  return elements_storage_.size();
}

template <typename Id, typename ES, typename... Fs>
Id ElementsContainer<Id, ES, Fs...>::Append() {
  Id result{GetCount()};
  elements_storage_.push_back({});
  features_storage_.push_back({});
  return result;
}

template <typename Id, typename ES, typename... Fs>
void ElementsContainer<Id, ES, Fs...>::Add(Id id) {
  std::ignore = GetOrInsert(elements_storage_, id);
  std::ignore = GetOrInsert(features_storage_, id);
}

template <typename Id, typename ES, typename... Fs>
void ElementsContainer<Id, ES, Fs...>::Initialize(size_t size) {
  elements_storage_.resize(size);
  features_storage_.resize(size);
}

template <typename Id, typename ES, typename... Fs>
template <typename F>
auto& ElementsContainer<Id, ES, Fs...>::GetFeatureStorage(Id id) {
  if constexpr (tuple_contains_v<std::tuple<Fs...>, F>) {
    return std::get<F>(features_storage_.at(id.value));
  } else {
    return elements_storage_.at(id.value).template GetFeatureStorage<F>();
  }
}

template <typename Id, typename ES, typename... Fs>
template <typename F>
const auto& ElementsContainer<Id, ES, Fs...>::GetFeatureStorage(Id id) const {
  if constexpr (tuple_contains_v<std::tuple<Fs...>, F>) {
    return std::get<F>(features_storage_.at(id.value));
  } else {
    return elements_storage_.at(id.value).template GetFeatureStorage<F>();
  }
}

template <typename Id, typename ES, typename... Fs>
template <typename F>
auto& ElementsContainer<Id, ES, Fs...>::GetFeatureExtraStorage() {
  if constexpr (tuple_contains_v<decltype(extra_features_storage_), F>) {
    return std::get<ExtraFeatureStorage<F>>(extra_features_storage_);
  } else {
    return std::get<ExtraFeatureStorage<F>>(elements_extra_features_storage_);
  }
}

template <typename Id, typename ES, typename... Fs>
template <typename F>
const auto& ElementsContainer<Id, ES, Fs...>::GetFeatureExtraStorage() const {
  if constexpr (tuple_contains_v<decltype(extra_features_storage_), F>) {
    return std::get<ExtraFeatureStorage<F>>(extra_features_storage_);
  } else {
    return std::get<ExtraFeatureStorage<F>>(elements_extra_features_storage_);
  }
}

#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename... Fs>
template <typename Feature>
inline constexpr bool ElementStorage<Fs...>::contains_element_feature =
    tuple_contains_v<std::tuple<Fs...>, Feature>;

template <typename... Fs>
template <typename F>
auto& ElementStorage<Fs...>::GetFeatureStorage() {
  return std::get<F>(features_storage_);
}

template <typename... Fs>
template <typename F>
const auto& ElementStorage<Fs...>::GetFeatureStorage() const {
  return std::get<F>(features_storage_);
}
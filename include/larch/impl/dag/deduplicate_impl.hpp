#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename CRTP, typename Feature>
auto& GetFeatureStorage(
    const FeatureConstView<Feature, CRTP, Deduplicate<Feature>>* feature) {
  const Feature* result = static_cast<const CRTP&>(*feature)
                              .template GetFeatureStorage<Deduplicate<Feature>>()
                              .feature_;
  Assert(result);
  return *result;
}

template <typename Feature, typename CRTP>
auto& FeatureMutableView<Deduplicate<Feature>, CRTP>::operator=(Feature&& feature) {
  auto& deduplicated = static_cast<CRTP&>(*this)
                           .template GetFeatureExtraStorage<Deduplicate<Feature>>()
                           .deduplicated_;
  const Feature* result =
      std::addressof(*deduplicated.insert(std::forward<Feature>(feature)).first);
  static_cast<CRTP&>(*this)
      .template GetFeatureStorage<Deduplicate<Feature>>()
      .feature_ = result;
  return *this;
}
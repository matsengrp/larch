#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Feature>
const Feature ExtraFeatureStorage<Deduplicate<Feature>>::empty_ = {};

template <typename CRTP, typename Feature>
auto& GetFeatureStorage(
    const FeatureConstView<Feature, CRTP, Deduplicate<Feature>>* feature) {
  const Feature* result = static_cast<const CRTP&>(*feature)
                              .template GetFeatureStorage<Deduplicate<Feature>>()
                              .feature_;
  // TODO Assert(result);
  if (not result) {
    return ExtraFeatureStorage<Deduplicate<Feature>>::empty_;
  }
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

template <typename Feature, typename CRTP>
const Feature* ExtraFeatureConstView<Deduplicate<Feature>, CRTP>::FindDeduplicated(
    const Feature& feature) const {
  using Id = std::conditional_t<
      CRTP::template contains_element_feature<NodeId, Deduplicate<Feature>>, NodeId,
      EdgeId>;
  auto& deduplicated = static_cast<const CRTP&>(*this)
                           .template GetFeatureExtraStorage<Id, Deduplicate<Feature>>()
                           .deduplicated_;
  auto result = deduplicated.find(feature);  // TODO Lock?
  if (result == deduplicated.end()) {
    return nullptr;
  }
  return std::addressof(*result);
}

template <typename Feature, typename CRTP>
std::pair<const Feature*, bool>
ExtraFeatureMutableView<Deduplicate<Feature>, CRTP>::AddDeduplicated(
    Feature&& feature) {
  using Id = std::conditional_t<
      CRTP::template contains_element_feature<NodeId, Deduplicate<Feature>>, NodeId,
      EdgeId>;
  auto& deduplicated = static_cast<CRTP&>(*this)
                           .template GetFeatureExtraStorage<Id, Deduplicate<Feature>>()
                           .deduplicated_;
  auto [iter, success] = deduplicated.insert(std::forward<Feature>(feature));
  return {std::addressof(*iter), success};
}
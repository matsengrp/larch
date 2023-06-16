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
// NOLINTNEXTLINE(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
auto& FeatureMutableView<Deduplicate<Feature>, CRTP>::operator=(
    Feature&& feature) const {
  auto& deduplicated = static_cast<const CRTP&>(*this)
                           .template GetFeatureExtraStorage<Deduplicate<Feature>>()
                           .deduplicated_;
  const Feature* result =
      std::addressof(*deduplicated.insert(std::forward<Feature>(feature)).first);
  static_cast<const CRTP&>(*this)
      .template GetFeatureStorage<Deduplicate<Feature>>()
      .feature_ = result;
  return *this;
}

template <typename Feature, typename CRTP>
// NOLINTNEXTLINE(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
auto& FeatureMutableView<Deduplicate<Feature>, CRTP>::operator=(
    const Feature* feature) const {
  static_cast<const CRTP&>(*this)
      .template GetFeatureStorage<Deduplicate<Feature>>()
      .feature_ = feature;
  return *this;
}

template <typename Feature, typename CRTP>
const Feature* ExtraFeatureConstView<Deduplicate<Feature>, CRTP>::FindDeduplicated(
    const Feature& feature) const {
  constexpr auto C =
      CRTP::template contains_element_feature<Component::Node, Deduplicate<Feature>>
          ? Component::Node
          : Component::Edge;
  auto& deduplicated = static_cast<const CRTP&>(*this)
                           .template GetFeatureExtraStorage<C, Deduplicate<Feature>>()
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
    const Feature& feature) const {
  constexpr auto C =
      CRTP::template contains_element_feature<Component::Node, Deduplicate<Feature>>
          ? Component::Node
          : Component::Edge;
  auto& deduplicated = static_cast<const CRTP&>(*this)
                           .template GetFeatureExtraStorage<C, Deduplicate<Feature>>()
                           .deduplicated_;
  auto existing = deduplicated.find(feature);
  if (existing != deduplicated.end()) {
    return {std::addressof(*existing), false};
  }
  auto [iter, success] = deduplicated.insert(feature.Copy());
  return {std::addressof(*iter), success};
}

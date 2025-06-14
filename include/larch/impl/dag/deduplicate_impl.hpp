#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Feature>
const Feature ExtraFeatureStorage<Deduplicate<Feature>>::empty_ = {};

template <typename CRTP, typename Feature>
auto GetFeatureStorage(
    const FeatureConstView<Feature, CRTP, Deduplicate<Feature>>* feature) {
  const Feature* result = static_cast<const CRTP&>(*feature)
                              .template GetFeatureStorage<Deduplicate<Feature>>()
                              .get()
                              .feature_;
  // TODO Assert(result);
  if (not result) {
    return std::cref(ExtraFeatureStorage<Deduplicate<Feature>>::empty_);
  }
  return std::cref(*result);
}

template <typename Feature, typename CRTP>
// NOLINTNEXTLINE(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
auto& FeatureMutableView<Deduplicate<Feature>, CRTP>::operator=(
    Feature&& feature) const {
  auto& deduplicated = static_cast<const CRTP&>(*this)
                           .template GetFeatureExtraStorage<Deduplicate<Feature>>()
                           .get()
                           .deduplicated_;
  const Feature* result =
      std::addressof(deduplicated.insert(std::forward<Feature>(feature)).first);
  static_cast<const CRTP&>(*this)
      .template GetFeatureStorage<Deduplicate<Feature>>()
      .get()
      .feature_ = result;
  return *this;
}

template <typename Feature, typename CRTP>
// NOLINTNEXTLINE(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
auto& FeatureMutableView<Deduplicate<Feature>, CRTP>::operator=(
    const Feature* feature) const {
  static_cast<const CRTP&>(*this)
      .template GetFeatureStorage<Deduplicate<Feature>>()
      .get()
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
                           .get()
                           .deduplicated_;
  auto result = deduplicated.find(feature);
  if (not result.has_value()) {
    return nullptr;
  } else {
    return std::addressof(result.value());
  }
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
                           .get()
                           .deduplicated_;

  auto result = deduplicated.insert(feature, [this](const auto& x) {
    return x.Copy(static_cast<const CRTP*>(this));
  });
  return {std::addressof(result.first), result.second};
}

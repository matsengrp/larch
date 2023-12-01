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
  const Feature* result = deduplicated.Write(
      [](auto& write, Feature&& feat) {
        return std::addressof(*write.insert(std::forward<Feature>(feat)).first);
      },
      std::forward<Feature>(feature));
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
  return deduplicated.Read(
      [](auto& read, const Feature& feat) {
        auto result = read.find(feat);
        if (result == read.end()) {
          return nullptr;
        }
        return std::addressof(*result);
      },
      feature);
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
  auto* existing = deduplicated.Read(
      [](auto& read, const Feature& feat) -> const Feature* {
        auto result = read.find(feat);
        if (result == read.end()) {
          return nullptr;
        }
        return std::addressof(*result);
      },
      feature);
  if (existing != nullptr) {
    return {existing, false};
  }
  return deduplicated.Write(
      [](auto& write, const Feature& feat) -> std::pair<const Feature*, bool> {
        auto [iter, success] = write.insert(feat.Copy());
        return {std::addressof(*iter), success};
      },
      feature);
}

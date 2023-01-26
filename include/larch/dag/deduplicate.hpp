#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename T>
using ConcurrentUnorderedSet =
    tbb::concurrent_unordered_set<T, std::hash<T>, std::equal_to<T>>;

/**
 * Used with any per-element feature to ensure that a single unique copy of
 * the undrlying feature is stored, and elements only store a pointer to it.
 * Stored pointers are stable.
 */
template <typename Feature>
struct Deduplicate {
  const Feature* feature_ = nullptr;
};

template <typename Feature>
struct ExtraFeatureStorage<Deduplicate<Feature>> {
  ExtraFeatureStorage() = default;
  MOVE_ONLY(ExtraFeatureStorage);

 private:
  template <typename CRTP, typename F>
  // NOLINTNEXTLINE(readability-redundant-declaration)
  friend auto& GetFeatureStorage(const FeatureConstView<F, CRTP, Deduplicate<F>>*);

  template <typename, typename, typename>
  friend struct FeatureMutableView;

  template <typename, typename>
  friend struct ExtraFeatureMutableView;
  ConcurrentUnorderedSet<Feature> deduplicated_;
  static const Feature empty_;
};

template <typename CRTP, typename Feature>
auto& GetFeatureStorage(
    const FeatureConstView<Feature, CRTP, Deduplicate<Feature>>* feature);

template <typename Feature, typename CRTP>
struct FeatureConstView<Deduplicate<Feature>, CRTP>
    : FeatureConstView<Feature, CRTP, Deduplicate<Feature>> {};

/**
 * Disables all modyfing functionality of the feature, and only allows to
 * fully replace the feature by the assignment operator.
 */
template <typename Feature, typename CRTP>
struct FeatureMutableView<Deduplicate<Feature>, CRTP>
    : FeatureConstView<Feature, CRTP, Deduplicate<Feature>> {
  auto& operator=(Feature&& feature);
  auto& operator=(const Feature* feature);
};

template <typename Feature, typename CRTP>
struct ExtraFeatureConstView<Deduplicate<Feature>, CRTP> {
  const Feature* FindDeduplicated(const Feature& feature) const;
};

template <typename Feature, typename CRTP>
struct ExtraFeatureMutableView<Deduplicate<Feature>, CRTP>
    : FeatureConstView<Feature, CRTP, Deduplicate<Feature>> {
  std::pair<const Feature*, bool> AddDeduplicated(Feature&& feature);
};
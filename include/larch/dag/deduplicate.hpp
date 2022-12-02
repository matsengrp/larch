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
  ConcurrentUnorderedSet<Feature> deduplicated_;
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
struct FeatureMutableView<Deduplicate<Feature>, CRTP> {
  auto& operator=(Feature&& feature);
};
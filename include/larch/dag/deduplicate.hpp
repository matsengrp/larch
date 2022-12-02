#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename T>
using ConcurrentUnorderedSet =
    tbb::concurrent_unordered_set<T, std::hash<T>, std::equal_to<T>>;

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

template <typename Feature, typename CRTP>
struct FeatureMutableView<Deduplicate<Feature>, CRTP> {
  auto& operator=(Feature&& feature);
};
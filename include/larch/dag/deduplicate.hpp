#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

#include <tbb/concurrent_unordered_set.h>

template <typename Feature>
class Deduplicate {
 public:
  inline static const Feature Empty = {};

  template <typename Id>
  class GlobalData {
   public:
    using DeduplicatedStorage =
        tbb::concurrent_unordered_set<Feature, std::hash<Feature>,
                                      std::equal_to<Feature>>;

    GlobalData(const GlobalData&) = delete;
    GlobalData(GlobalData&&) = default;
    GlobalData& operator=(const GlobalData&) = delete;
    GlobalData& operator=(GlobalData&&) = default;

    GlobalData();

    const Feature& Get(const Deduplicate<Feature>& self, Id id) const;
    void Set(Deduplicate<Feature>& self, Id id, Feature&& feature);

    const Feature* ContainsUnique(const Feature& feature) const;
    const Feature* AddUnique(Feature&& feature);

   private:
    DeduplicatedStorage deduplicated_storage_;
  };

 private:
  const Feature* feature_ = &Empty;
};

template <typename Feature, typename View>
struct FeatureTraits<Deduplicate<Feature>, View> {
  using Reader = FeatureReader<Feature, View>;
  using Writer = FeatureWriter<Feature, View>;
  template <typename Id>
  using GlobalData = typename Deduplicate<Feature>::GlobalData<Id>;
};
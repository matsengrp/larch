#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

#include <tbb/concurrent_unordered_map.h>

template <typename Id, typename Feature>
class Deduplicate {
 public:
  using Map =
      tbb::concurrent_unordered_map<Id, Feature, std::hash<Id>, std::equal_to<Id>>;

  class GlobalData {
   public:
   private:
    Map map_;
  };

 private:
  Feature* feature_;
};

template <typename Id, typename Feature, typename View>
class FeatureReader<Deduplicate<Id, Feature>, View> {
 public:
};

template <typename Id, typename Feature, typename View>
class FeatureWriter<Deduplicate<Id, Feature>, View>
    : public FeatureReader<Deduplicate<Id, Feature>, View> {
 public:
};
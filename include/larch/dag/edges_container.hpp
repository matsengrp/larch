#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Storage, typename... Features>
class DefaultEdgesContainer {
 public:
  using StorageType = Storage;
  using FeaturesType = std::tuple<Features...>;

  template <typename DAG>
  class ViewBase
      : public std::conditional_t<DAG::is_mutable,
                                  FeatureWriter<Features, EdgeView<DAG>>,
                                  FeatureReader<Features, EdgeView<DAG>>>... {};

 private:
  DAG_VIEW_FRIENDS;
  std::vector<Storage> edges_;
  [[no_unique_address]] std::tuple<Features...> features_;
};
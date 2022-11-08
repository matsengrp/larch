#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Storage, typename... Features>
class DefaultNodesContainer {
 public:
  using StorageType = Storage;
  using FeaturesType = std::tuple<Features...>;

  template <typename DAG>
  class ViewBase
      : public std::conditional_t<DAG::is_mutable,
                                  FeatureWriter<Features, NodeView<DAG>>,
                                  FeatureReader<Features, NodeView<DAG>>>... {};

 private:
  DAG_VIEW_FRIENDS;
  std::vector<Storage> nodes_;
  [[no_unique_address]] std::tuple<Features...> features_;
};
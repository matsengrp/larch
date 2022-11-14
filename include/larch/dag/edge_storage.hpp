#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename... Features>
class DefaultEdgeStorage {
 public:
  template <typename DAG>
  using ViewType = EdgeView<DAG>;

  template <typename DAG>
  class ViewBase
      : public std::conditional_t<
            DAG::is_mutable, typename FeatureTraits<Features, ViewType<DAG>>::Writer,
            typename FeatureTraits<Features, ViewType<DAG>>::Reader>... {};

 private:
  DAG_FEATURE_FRIENDS;
  DAG_VIEW_FRIENDS;
  NodeId parent_;
  NodeId child_;
  CladeIdx clade_;
  [[no_unique_address]] std::tuple<Features...> features_;
};
#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename... Features>
class DefaultEdgeStorage {
 public:
  template <typename DAG>
  using ViewType = EdgeView<DAG, Features...>;

 private:
  DAG_FEATURE_FRIENDS;
  DAG_VIEW_FRIENDS;
  NodeId parent_;
  NodeId child_;
  CladeIdx clade_;
  [[no_unique_address]] std::tuple<Features...> features_;
};
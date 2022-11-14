#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename NodesContainer, typename EdgesContainer, typename... Features>
class DefaultDAGStorage {
 public:
  using NodesContainerType = NodesContainer;
  using EdgesContainerType = EdgesContainer;

  template <typename Storage>
  class DAGViewBase
      : public std::conditional_t<std::is_const_v<Storage>,
                                  FeatureReader<Features, DAGView<Storage>>,
                                  FeatureWriter<Features, DAGView<Storage>>>... {};

  DefaultDAGStorage() = default;
  MOVE_ONLY(DefaultDAGStorage);

  auto View();
  auto View() const;

 private:
  DAG_FEATURE_FRIENDS;
  DAG_VIEW_FRIENDS;
  NodesContainer nodes_;
  EdgesContainer edges_;
  NodeId root_ = {NoId};
  std::vector<NodeId> leafs_;
  [[no_unique_address]] std::tuple<Features...> features_;
};
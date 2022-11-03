#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename NodesContainer, typename EdgesContainer, typename... Features>
class DefaultDAGStorage {
 public:
  using NodesContainerType = NodesContainer;
  using EdgesContainerType = EdgesContainer;
  
  DefaultDAGStorage() = default;
  MOVE_ONLY(DefaultDAGStorage);

  auto View();

 private:
  DAG_FEATURE_FRIENDS;
  DAG_VIEW_FRIENDS;
  NodesContainer nodes_;
  EdgesContainer edges_;
  [[no_unique_address]] std::tuple<Features...> features_;
};
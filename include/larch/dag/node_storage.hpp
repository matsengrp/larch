#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename... Features>
class DefaultNodeStorage {
 public:
  template <typename DAG>
  using ViewType = NodeView<DAG, Features...>;

  void ClearConnections();
  void AddEdge(CladeIdx clade, EdgeId id, bool this_node_is_parent);

 private:
  DAG_FEATURE_FRIENDS;
  DAG_VIEW_FRIENDS;
  std::vector<EdgeId> parents_;
  std::vector<std::vector<EdgeId>> clades_;
  [[no_unique_address]] std::tuple<Features...> features_;
};
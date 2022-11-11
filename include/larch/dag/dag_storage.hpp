#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename NodesContainer, typename EdgesContainer, typename... Features>
class DefaultDAGStorage {
 public:
  using NodesContainerType = NodesContainer;
  using EdgesContainerType = EdgesContainer;

  using PerDAGFeatures = filter_t<feature_is_per_dag, std::tuple, Features...>;
  using PerNodeFeatures = filter_t<feature_is_per_node, std::tuple, Features...>;
  using PerEdgeFeatures = filter_t<feature_is_per_edge, std::tuple, Features...>;

  template <typename Feature>
  static inline constexpr bool contains_node_feature =
      tuple_contians_v<PerNodeFeatures, Feature>;

  template <typename Feature>
  static inline constexpr bool contains_edge_feature =
      tuple_contians_v<PerEdgeFeatures, Feature>;

  template <typename Storage, typename... Fs>
  class DAGViewBase
      : public std::conditional_t<std::is_const_v<Storage>,
                                  FeatureReader<Fs, DAGView<Storage>>,
                                  FeatureWriter<Fs, DAGView<Storage>>>... {};

  template <typename Storage, typename... Fs>
  class NodeViewBase
      : public std::conditional_t<std::is_const_v<Storage>,
                                  FeatureReader<Fs, NodeView<Storage>>,
                                  FeatureWriter<Fs, NodeView<Storage>>>... {};

  template <typename Storage, typename... Fs>
  class EdgeViewBase
      : public std::conditional_t<std::is_const_v<Storage>,
                                  FeatureReader<Fs, EdgeView<Storage>>,
                                  FeatureWriter<Fs, EdgeView<Storage>>>... {};

  template <typename Storage, typename... Fs>
  static constexpr auto GetDAGViewBase(std::tuple<Fs...>) {
    return DAGViewBase<Storage, Fs...>{};
  }

  template <typename Storage, typename... Fs>
  static constexpr auto GetNodeViewBase(std::tuple<Fs...>) {
    return NodeViewBase<Storage, Fs...>{};
  }

  template <typename Storage, typename... Fs>
  static constexpr auto GetEdgeViewBase(std::tuple<Fs...>) {
    return EdgeViewBase<Storage, Fs...>{};
  }

  template <typename Storage>
  using DAGViewBaseType = decltype(GetDAGViewBase<Storage>(PerDAGFeatures{}));

  template <typename Storage>
  using NodeViewBaseType = decltype(GetNodeViewBase<Storage>(PerNodeFeatures{}));

  template <typename Storage>
  using EdgeViewBaseType = decltype(GetEdgeViewBase<Storage>(PerEdgeFeatures{}));

  DefaultDAGStorage() = default;
  MOVE_ONLY(DefaultDAGStorage);

  auto View();

 private:
  DAG_FEATURE_FRIENDS;
  DAG_VIEW_FRIENDS;
  NodesContainer nodes_;
  EdgesContainer edges_;
  NodeId root_ = {NoId};
  std::vector<NodeId> leafs_;
  [[no_unique_address]] std::tuple<Features...> features_;
};
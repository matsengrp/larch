#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Helper types for use with ExtendDAGStorage to add more features
 * to an existing DAG. Example use:
 *
 *     ExtendDAGStorage my_extended_dag{existing_dag_view,
 *         Extend::Nodes<Deduplicate<MyNewNodeFeature>>{},
 *         Extend::Edges<MyNewEdgeFeature1, MyNewEdgeFeature2>{}};
 *
 */
namespace Extend {

template <typename... Fs>
struct Nodes {
  using FeatureTypes = std::tuple<Fs...>;
  using Storage = std::vector<std::tuple<Fs...>>;
  using ExtraStorage = std::tuple<ExtraFeatureStorage<Fs>...>;
  template <typename Self, typename CRTP>
  struct ConstView : FeatureConstView<Fs, CRTP>... {};
  template <typename Self, typename CRTP>
  struct MutableView : ConstView<Self, CRTP>, FeatureMutableView<Fs, CRTP>... {
    using FeatureMutableView<Fs, CRTP>::operator=...;
  };

  template <typename Self, typename CRTP>
  struct ExtraConstView : ExtraFeatureConstView<Fs, CRTP>... {};
  template <typename Self, typename CRTP>
  struct ExtraMutableView : ExtraConstView<Self, CRTP>,
                            ExtraFeatureMutableView<Fs, CRTP>... {};
  template <typename Feature>
  static const bool contains_element_feature =
      tuple_contains_v<std::tuple<Fs...>, Feature>;
};

template <typename... Fs>
struct Edges {
  using FeatureTypes = std::tuple<Fs...>;
  using Storage = std::vector<std::tuple<Fs...>>;
  using ExtraStorage = std::tuple<ExtraFeatureStorage<Fs>...>;
  template <typename Self, typename CRTP>
  struct ConstView : FeatureConstView<Fs, CRTP>... {};
  template <typename Self, typename CRTP>
  struct MutableView : ConstView<Self, CRTP>, FeatureMutableView<Fs, CRTP>... {
    using FeatureMutableView<Fs, CRTP>::operator=...;
  };

  template <typename Self, typename CRTP>
  struct ExtraConstView : ExtraFeatureConstView<Fs, CRTP>... {};
  template <typename Self, typename CRTP>
  struct ExtraMutableView : ExtraConstView<Self, CRTP>,
                            ExtraFeatureMutableView<Fs, CRTP>... {};

  template <typename Feature>
  static const bool contains_element_feature =
      tuple_contains_v<std::tuple<Fs...>, Feature>;
};

template <typename... Fs>
struct DAG {
  using FeatureTypes = std::tuple<Fs...>;
  using Storage = std::tuple<Fs...>;
  template <typename Self, typename CRTP>
  struct ConstView : FeatureConstView<Fs, CRTP>...,
                     ExtraFeatureConstView<Fs, CRTP>... {};
  template <typename Self, typename CRTP>
  struct MutableView : ConstView<Self, CRTP>,
                       FeatureMutableView<Fs, CRTP>...,
                       ExtraFeatureMutableView<Fs, CRTP>... {};
};

template <typename...>
struct Empty {
  using FeatureTypes = std::tuple<>;
  using Storage = std::tuple<>;
  template <typename, typename>
  struct ConstView {};
  template <typename, typename>
  struct MutableView {};
  struct ExtraConstView {};
  template <typename, typename>
  struct ExtraMutableView {};
  template <typename>
  static const bool contains_element_feature = false;
};
}  // namespace Extend

/**
 * Adds new features to an existing DAG. See `namespace Extend` for more info.
 */
template <typename Target, typename Arg0 = Extend::Empty<>,
          typename Arg1 = Extend::Empty<>, typename Arg2 = Extend::Empty<>>
struct ExtendDAGStorage {
 public:
  using Self = ExtendDAGStorage<Target, Arg0, Arg1, Arg2>;
  using TargetView = decltype(ViewOf(std::declval<Target>()));
  using OnNodes = select_argument_t<Extend::Nodes, Arg0, Arg1, Arg2>;
  using OnEdges = select_argument_t<Extend::Edges, Arg0, Arg1, Arg2>;
  using OnDAG = select_argument_t<Extend::DAG, Arg0, Arg1, Arg2>;

  using FeatureTypes = typename TargetView::StorageType::FeatureTypes;
  using AllNodeFeatures =
      decltype(std::tuple_cat(typename OnNodes::FeatureTypes{},
                              typename TargetView::StorageType::AllNodeFeatures{}));
  using AllEdgeFeatures =
      decltype(std::tuple_cat(typename OnEdges::FeatureTypes{},
                              typename TargetView::StorageType::AllEdgeFeatures{}));

  template <typename Id, typename CRTP>
  struct ConstElementViewBase;

  template <typename CRTP>
  struct ConstElementViewBase<NodeId, CRTP>
      : TargetView::StorageType::template ConstElementViewBase<NodeId, CRTP>,
        OnNodes::template ConstView<Self, CRTP> {};

  template <typename CRTP>
  struct ConstElementViewBase<EdgeId, CRTP>
      : TargetView::StorageType::template ConstElementViewBase<EdgeId, CRTP>,
        OnEdges::template ConstView<Self, CRTP> {};

  template <typename Id, typename CRTP>
  struct MutableElementViewBase;

  template <typename CRTP>
  struct MutableElementViewBase<NodeId, CRTP>
      : TargetView::StorageType::template MutableElementViewBase<NodeId, CRTP>,
        OnNodes::template MutableView<Self, CRTP> {
    using TargetView::StorageType::template MutableElementViewBase<NodeId,
                                                                   CRTP>::operator=;
    using OnNodes::template MutableView<Self, CRTP>::operator=;
  };

  template <typename CRTP>
  struct MutableElementViewBase<EdgeId, CRTP>
      : TargetView::StorageType::template MutableElementViewBase<EdgeId, CRTP>,
        OnEdges::template MutableView<Self, CRTP> {
    using TargetView::StorageType::template MutableElementViewBase<EdgeId,
                                                                   CRTP>::operator=;
    using OnEdges::template MutableView<Self, CRTP>::operator=;
  };

  template <typename Id, typename Feature>
  static const bool contains_element_feature;

  template <typename CRTP>
  struct ConstDAGViewBase : TargetView::StorageType::template ConstDAGViewBase<CRTP>,
                            OnDAG::template ConstView<Self, CRTP>,
                            OnNodes::template ExtraConstView<Self, CRTP>,
                            OnEdges::template ExtraConstView<Self, CRTP> {};
  template <typename CRTP>
  struct MutableDAGViewBase
      : TargetView::StorageType::template MutableDAGViewBase<CRTP>,
        OnDAG::template MutableView<Self, CRTP>,
        OnNodes::template ExtraMutableView<Self, CRTP>,
        OnEdges::template ExtraMutableView<Self, CRTP> {};

  MOVE_ONLY(ExtendDAGStorage);

  ExtendDAGStorage();
  explicit ExtendDAGStorage(Target&& target);

  auto View();
  auto View() const;

  NodeId AppendNode();
  EdgeId AppendEdge();

  void AddNode(NodeId id);
  void AddEdge(EdgeId id);

  size_t GetNodesCount() const;
  size_t GetEdgesCount() const;

  auto GetNodes() const;
  auto GetEdges() const;

  void InitializeNodes(size_t size);
  void InitializeEdges(size_t size);

  template <typename F>
  auto& GetFeatureStorage();

  template <typename F>
  const auto& GetFeatureStorage() const;

  template <typename F>
  auto& GetFeatureStorage(NodeId id);

  template <typename F>
  const auto& GetFeatureStorage(NodeId id) const;

  template <typename F>
  auto& GetFeatureStorage(EdgeId id);

  template <typename F>
  const auto& GetFeatureStorage(EdgeId id) const;

  template <typename Id, typename F>
  auto& GetFeatureExtraStorage();

  template <typename Id, typename F>
  const auto& GetFeatureExtraStorage() const;

 private:
  auto GetTarget();
  auto GetTarget() const;

  std::decay_t<Target> target_ = {};
  typename OnNodes::Storage additional_node_features_storage_;
  typename OnEdges::Storage additional_edge_features_storage_;
  typename OnDAG::Storage additional_dag_features_storage_;
  typename OnNodes::ExtraStorage additional_node_extra_features_storage_;
  typename OnEdges::ExtraStorage additional_edge_extra_features_storage_;
};

template <typename Target, typename Arg0 = Extend::Empty<>,
          typename Arg1 = Extend::Empty<>, typename Arg2 = Extend::Empty<>>
auto ExtendStorage(Target&& target, Arg0 = Extend::Empty<>{}, Arg1 = Extend::Empty<>{},
                   Arg2 = Extend::Empty<>{});

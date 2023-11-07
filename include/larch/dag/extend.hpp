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
  struct MutableView : FeatureMutableView<Fs, CRTP>... {
    using FeatureMutableView<Fs, CRTP>::operator=...;
  };

  template <typename Self, typename CRTP>
  struct ExtraConstView : ExtraFeatureConstView<Fs, CRTP>... {};
  template <typename Self, typename CRTP>
  struct ExtraMutableView : ExtraFeatureMutableView<Fs, CRTP>... {};
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
  struct MutableView : FeatureMutableView<Fs, CRTP>... {
    using FeatureMutableView<Fs, CRTP>::operator=...;
  };

  template <typename Self, typename CRTP>
  struct ExtraConstView : ExtraFeatureConstView<Fs, CRTP>... {};
  template <typename Self, typename CRTP>
  struct ExtraMutableView : ExtraFeatureMutableView<Fs, CRTP>... {};

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
  struct MutableView : FeatureMutableView<Fs, CRTP>...,
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
template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
struct ExtendDAGStorage {
 public:
  constexpr static const Component component = Component::DAG;
  constexpr static const Role role = Role::Storage;

  static_assert(not std::is_reference_v<Target>);
  static_assert(Target::component == Component::DAG);

  static_assert(IsNameCorrect<ShortName, ExtendDAGStorage>::value);

  using Self =
      std::conditional_t<std::is_same_v<ShortName, void>, ExtendDAGStorage, ShortName>;

  using TargetView = typename ViewTypeOf<Target>::type;
  using OnNodes = select_argument_t<Extend::Nodes, Arg0, Arg1, Arg2>;
  using OnEdges = select_argument_t<Extend::Edges, Arg0, Arg1, Arg2>;
  using OnDAG = select_argument_t<Extend::DAG, Arg0, Arg1, Arg2>;

  struct ExtraStorageType {
    using FeatureTypes = decltype(std::tuple_cat(
        typename TargetView::StorageType::ExtraStorageType::FeatureTypes{},
        typename OnDAG::FeatureTypes{}));

    template <template <typename, typename> typename T, typename CRTP>
    struct Base : TargetView::StorageType::ExtraStorageType::template Base<T, CRTP> {};

    template <typename CRTP>
    struct Base<FeatureConstView, CRTP>
        : TargetView::StorageType::template ConstDAGViewBase<CRTP>,
          OnDAG::template ConstView<Self, CRTP> {};

    template <typename CRTP>
    struct Base<FeatureMutableView, CRTP>
        : TargetView::StorageType::template MutableDAGViewBase<CRTP>,
          OnDAG::template MutableView<Self, CRTP> {};
  };

  using FeatureTypes = typename TargetView::StorageType::FeatureTypes;
  template <Component C>
  using Container = typename TargetView::StorageType::template Container<C>;
  using AllNodeFeatures =
      decltype(std::tuple_cat(typename OnNodes::FeatureTypes{},
                              typename TargetView::StorageType::AllNodeFeatures{}));
  using AllEdgeFeatures =
      decltype(std::tuple_cat(typename OnEdges::FeatureTypes{},
                              typename TargetView::StorageType::AllEdgeFeatures{}));

  template <Component C, typename CRTP>
  struct ConstElementViewBase;

  template <typename CRTP>
  struct ConstElementViewBase<Component::Node, CRTP>
      : TargetView::StorageType::template ConstElementViewBase<Component::Node, CRTP>,
        OnNodes::template ConstView<Self, CRTP> {};

  template <typename CRTP>
  struct ConstElementViewBase<Component::Edge, CRTP>
      : TargetView::StorageType::template ConstElementViewBase<Component::Edge, CRTP>,
        OnEdges::template ConstView<Self, CRTP> {};

  template <Component C, typename CRTP>
  struct MutableElementViewBase;

  template <typename CRTP>
  struct MutableElementViewBase<Component::Node, CRTP>
      : TargetView::StorageType::template MutableElementViewBase<Component::Node, CRTP>,
        OnNodes::template MutableView<Self, CRTP> {
    using TargetView::StorageType::template MutableElementViewBase<Component::Node,
                                                                   CRTP>::operator=;
    using OnNodes::template MutableView<Self, CRTP>::operator=;
  };

  template <typename CRTP>
  struct MutableElementViewBase<Component::Edge, CRTP>
      : TargetView::StorageType::template MutableElementViewBase<Component::Edge, CRTP>,
        OnEdges::template MutableView<Self, CRTP> {
    using TargetView::StorageType::template MutableElementViewBase<Component::Edge,
                                                                   CRTP>::operator=;
    using OnEdges::template MutableView<Self, CRTP>::operator=;
  };

  template <Component C, typename Feature>
  static const bool contains_element_feature;

  template <typename CRTP>
  struct ConstDAGViewBase : ExtraStorageType::template Base<FeatureConstView, CRTP>,
                            OnNodes::template ExtraConstView<Self, CRTP>,
                            OnEdges::template ExtraConstView<Self, CRTP> {};
  template <typename CRTP>
  struct MutableDAGViewBase : ExtraStorageType::template Base<FeatureMutableView, CRTP>,
                              OnNodes::template ExtraMutableView<Self, CRTP>,
                              OnEdges::template ExtraMutableView<Self, CRTP> {
    template <typename F>
    constexpr auto& AsFeature() const noexcept {
      return static_cast<const ExtraFeatureMutableView<F, CRTP>&>(*this);
    }
  };

  MOVE_ONLY(ExtendDAGStorage);

  static Self Consume(Target&& target) {
    static_assert(Target::role == Role::Storage);
    return Self{std::move(target)};
  }

  static Self FromView(const Target& target) {
    static_assert(Target::role == Role::View);
    return Self{Target{target}};
  }

  static Self EmptyDefault() {
    static_assert(Target::role == Role::Storage);
    return Self{Target{}};
  }

  DAGView<Self> View();
  DAGView<const Self> View() const;

  NodeId AppendNode();
  EdgeId AppendEdge();

  void AddNode(NodeId id);
  void AddEdge(EdgeId id);

  size_t GetNodesCount() const;
  size_t GetEdgesCount() const;

  template <Component C>
  Id<C> GetNextAvailableId() const {
    return GetTarget().template GetNextAvailableId<C>();
  }

  auto GetNodes() const;
  auto GetEdges() const;

  void InitializeNodes(size_t size);
  void InitializeEdges(size_t size);

  void ClearNodes();
  void ClearEdges();

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

  template <Component C, typename F>
  auto& GetFeatureExtraStorage();

  template <Component C, typename F>
  const auto& GetFeatureExtraStorage() const;

  template <Component C>
  auto GetContainer() -> typename TargetView::StorageType::template Container<C>&;

  template <Component C>
  auto GetContainer() const -> const
      typename TargetView::StorageType::template Container<C>&;

  auto& GetTargetStorage() { return *this; }
  auto& GetTargetStorage() const { return *this; }

 private:
  friend ShortName;
  explicit ExtendDAGStorage(Target&& target);

  auto GetTarget();
  auto GetTarget() const;

  Target target_;
  typename OnNodes::Storage additional_node_features_storage_;
  typename OnEdges::Storage additional_edge_features_storage_;
  typename OnDAG::Storage additional_dag_features_storage_;
  typename OnNodes::ExtraStorage additional_node_extra_features_storage_;
  typename OnEdges::ExtraStorage additional_edge_extra_features_storage_;
};

template <typename ShortName, typename Target, typename Arg0 = Extend::Empty<>,
          typename Arg1 = Extend::Empty<>, typename Arg2 = Extend::Empty<>>
using ExtendStorageType = ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>;

template <typename ShortName, typename Target, typename Arg0 = Extend::Empty<>,
          typename Arg1 = Extend::Empty<>, typename Arg2 = Extend::Empty<>,
          typename = std::enable_if_t<Target::role == Role::Storage>>
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2> AddExtend(Target&& target) {
  return ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::Consume(
      std::move(target));
}

template <typename ShortName, typename Target, typename Arg0 = Extend::Empty<>,
          typename Arg1 = Extend::Empty<>, typename Arg2 = Extend::Empty<>,
          typename = std::enable_if_t<Target::role == Role::View>>
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2> AddExtend(const Target& target) {
  return ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::FromView(target);
}

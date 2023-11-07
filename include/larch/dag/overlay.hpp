#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

struct Overlay {};

template <typename CRTP, typename Tag>
struct FeatureConstView<Overlay, CRTP, Tag> {
  template <typename F>
  bool IsOverlaid() const;
  bool IsAppended() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<Overlay, CRTP, Tag> {
  template <typename F>
  auto SetOverlay() const;
};

struct OverlayDAG {};  // TODO make it extra feature

template <typename CRTP, typename Tag>
struct FeatureConstView<OverlayDAG, CRTP, Tag> {
  auto GetOriginal() const;
  bool HaveOverlays() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<OverlayDAG, CRTP, Tag> {};

namespace {
template <typename>
struct ToOverlayStorage;

template <typename... Features>
struct ToOverlayStorage<std::tuple<Features...>> {
  template <typename Id>
  using type = std::tuple<std::unordered_map<Id, Features>...>;
};

}  // namespace

template <typename ShortName, typename Target>
struct OverlayDAGStorage {
  constexpr static const Component component = Component::DAG;
  constexpr static const Role role = Role::Storage;

  static_assert(not std::is_reference_v<Target>);
  static_assert(Target::component == Component::DAG);

  static_assert(IsNameCorrect<ShortName, OverlayDAGStorage>::value);

  using Self =
      std::conditional_t<std::is_same_v<ShortName, void>, OverlayDAGStorage, ShortName>;

  using TargetView = typename ViewTypeOf<Target>::type;

  struct ExtraStorageType : TargetView::StorageType::ExtraStorageType {
    MOVE_ONLY(ExtraStorageType);

    using FeatureTypes = decltype(std::tuple_cat(
        typename TargetView::StorageType::ExtraStorageType::FeatureTypes{},
        std::tuple<OverlayDAG>{}));

    template <template <typename, typename> typename T, typename CRTP>
    struct Base : TargetView::StorageType::ExtraStorageType::template Base<T, CRTP> {};

    template <typename CRTP>
    struct Base<FeatureConstView, CRTP>
        : TargetView::StorageType::template ConstDAGViewBase<CRTP>,
          FeatureConstView<OverlayDAG, CRTP> {};

    template <typename CRTP>
    struct Base<FeatureMutableView, CRTP>
        : TargetView::StorageType::template MutableDAGViewBase<CRTP>,
          FeatureMutableView<OverlayDAG, CRTP> {};
  };

  using FeatureTypes = typename ExtraStorageType::FeatureTypes;
  template <Component C>
  using Container = typename TargetView::StorageType::template Container<C>;
  using AllNodeFeatures = typename TargetView::StorageType::AllNodeFeatures;
  using AllEdgeFeatures = typename TargetView::StorageType::AllEdgeFeatures;

  template <Component C, typename CRTP>
  struct ConstElementViewBase
      : FeatureConstView<Overlay, CRTP>,
        TargetView::StorageType::template ConstElementViewBase<C, CRTP> {};

  template <Component C, typename CRTP>
  struct MutableElementViewBase
      : FeatureMutableView<Overlay, CRTP>,
        TargetView::StorageType::template MutableElementViewBase<C, CRTP> {
    using TargetView::StorageType::template MutableElementViewBase<C, CRTP>::operator=;
  };

  template <Component C, typename Feature>
  static const bool contains_element_feature;

  template <typename CRTP>
  struct ConstDAGViewBase : ExtraStorageType::template Base<FeatureConstView, CRTP> {};

  template <typename CRTP>
  struct MutableDAGViewBase
      : ExtraStorageType::template Base<FeatureMutableView, CRTP> {};

  MOVE_ONLY(OverlayDAGStorage);

  static Self Consume(Target&& target) {
    static_assert(Target::role == Role::Storage);
    return Self{std::move(target)};
  }

  static Self FromView(const Target& target) {
    static_assert(Target::role == Role::View);
    return Self{Target{target}};
  }

  auto View();
  auto View() const;

  NodeId AppendNode();
  EdgeId AppendEdge();

  void AddNode(NodeId id);
  void AddEdge(EdgeId id);

  size_t GetNodesCount() const;
  size_t GetEdgesCount() const;

  template <Component C>
  Id<C> GetNextAvailableId() const {
    if constexpr (C == Component::Node) {
      return {GetTarget().template GetNextAvailableId<C>().value +
              added_node_storage_.size()};
    } else {
      return {GetTarget().template GetNextAvailableId<C>().value +
              added_edge_storage_.size()};
    }
  }

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

  template <Component C, typename F>
  auto& GetFeatureExtraStorage();

  template <Component C, typename F>
  const auto& GetFeatureExtraStorage() const;

  auto& GetTargetStorage() { return *this; }
  auto& GetTargetStorage() const { return *this; }

 private:
  friend ShortName;
  explicit OverlayDAGStorage(Target&& target);

  auto GetTarget();
  auto GetTarget() const;

  template <typename F, typename OverlayStorageType>
  static auto GetFeatureStorageImpl(OverlayStorageType& self, NodeId id)
      -> std::conditional_t<not std::is_const_v<OverlayStorageType> and
                                OverlayStorageType::TargetView::is_mutable,
                            F&, const F&>;

  template <typename F, typename OverlayStorageType>
  static auto GetFeatureStorageImpl(OverlayStorageType& self, EdgeId id)
      -> std::conditional_t<not std::is_const_v<OverlayStorageType> and
                                OverlayStorageType::TargetView::is_mutable,
                            F&, const F&>;

  template <typename, typename, typename>
  friend struct FeatureConstView;

  template <typename, typename, typename>
  friend struct FeatureMutableView;

  Target target_;
  typename ToOverlayStorage<AllNodeFeatures>::template type<NodeId>
      replaced_node_storage_;
  typename ToOverlayStorage<AllEdgeFeatures>::template type<EdgeId>
      replaced_edge_storage_;
  std::vector<AllNodeFeatures> added_node_storage_;
  std::vector<AllEdgeFeatures> added_edge_storage_;
};

template <typename ShortName, typename DAG>
using OverlayStorageType = OverlayDAGStorage<ShortName, DAG>;

template <typename ShortName, typename DAG,
          typename = std::enable_if_t<DAG::role == Role::Storage>>
OverlayDAGStorage<ShortName, DAG> AddOverlay(DAG&& dag) {
  return OverlayDAGStorage<ShortName, DAG>::Consume(std::move(dag));
}

template <typename ShortName, typename DAG,
          typename = std::enable_if_t<DAG::role == Role::View>>
typename OverlayDAGStorage<ShortName, DAG>::Self AddOverlay(const DAG& dag) {
  return OverlayDAGStorage<ShortName, DAG>::FromView(dag);
}

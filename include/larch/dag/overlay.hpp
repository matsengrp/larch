#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

struct Overlay {};

template <typename CRTP, typename Tag>
struct FeatureConstView<Overlay, CRTP, Tag> {
  bool IsOverlaid() const;
  bool IsAppended() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<Overlay, CRTP, Tag> {
  void Overlay();
};

struct OverlayDAG {};  // TODO make it extra feature

template <typename CRTP, typename Tag>
struct FeatureConstView<OverlayDAG, CRTP, Tag> {
  auto GetOriginal() const;
  bool HaveOverlays() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<OverlayDAG, CRTP, Tag>
    : FeatureConstView<OverlayDAG, CRTP, Tag> {};

template <typename Target>
struct OverlayDAGStorage {
 public:
  using TargetView = decltype(ViewOf(std::declval<Target>()));

  using FeatureTypes = typename TargetView::StorageType::FeatureTypes;
  using AllNodeFeatures = typename TargetView::StorageType::AllNodeFeatures;
  using AllEdgeFeatures = typename TargetView::StorageType::AllEdgeFeatures;

  template <typename Id, typename CRTP>
  struct ConstElementViewBase
      : FeatureConstView<Overlay, CRTP>,
        TargetView::StorageType::template ConstElementViewBase<Id, CRTP> {};

  template <typename Id, typename CRTP>
  struct MutableElementViewBase
      : FeatureConstView<Overlay, CRTP>,
        FeatureMutableView<Overlay, CRTP>,
        TargetView::StorageType::template MutableElementViewBase<Id, CRTP> {
    using TargetView::StorageType::template MutableElementViewBase<Id, CRTP>::operator=;
  };

  template <typename Id, typename Feature>
  static const bool contains_element_feature;

  template <typename CRTP>
  struct ConstDAGViewBase : FeatureConstView<OverlayDAG, CRTP>,
                            TargetView::StorageType::template ConstDAGViewBase<CRTP> {};

  template <typename CRTP>
  struct MutableDAGViewBase
      : FeatureMutableView<OverlayDAG, CRTP>,
        TargetView::StorageType::template MutableDAGViewBase<CRTP> {};

  MOVE_ONLY(OverlayDAGStorage);

  explicit OverlayDAGStorage(Target&& target);

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

  // TODO write const/non-const helpers for the rest of feature getters
  template <typename F, typename OverlayStorageType>
  static auto GetFeatureStorageImpl(OverlayStorageType& self, NodeId id)
      -> std::conditional_t<OverlayStorageType::TargetView::is_mutable, F&, const F&>;

  template <typename, typename, typename>
  friend struct FeatureConstView;

  template <typename, typename, typename>
  friend struct FeatureMutableView;

  std::decay_t<Target> target_ = {};
  std::unordered_map<NodeId, AllNodeFeatures> replaced_node_storage_;
  std::unordered_map<EdgeId, AllEdgeFeatures> replaced_edge_storage_;
  std::vector<AllNodeFeatures> added_node_storage_;
  std::vector<AllEdgeFeatures> added_edge_storage_;
};

template <typename DAG>
auto AddOverlay(DAG&& dag) {
  return OverlayDAGStorage<DAG>{std::forward<DAG>(dag)};
}

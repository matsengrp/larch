#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Target>
struct OverlayDAGStorage {
 public:
  using TargetView = decltype(ViewOf(std::declval<Target>()));

  using FeatureTypes = typename TargetView::StorageType::FeatureTypes;
  using AllNodeFeatures = typename TargetView::StorageType::AllNodeFeatures;
  using AllEdgeFeatures = typename TargetView::StorageType::AllEdgeFeatures;

  template <typename Id, typename CRTP>
  struct ConstElementViewBase
      : TargetView::StorageType::template ConstElementViewBase<Id, CRTP> {};

  template <typename Id, typename CRTP>
  struct MutableElementViewBase
      : TargetView::StorageType::template MutableElementViewBase<Id, CRTP> {
    using TargetView::StorageType::template MutableElementViewBase<Id, CRTP>::operator=;
  };

  template <typename Id, typename Feature>
  static const bool contains_element_feature;

  template <typename CRTP>
  struct ConstDAGViewBase : TargetView::StorageType::template ConstDAGViewBase<CRTP> {};

  template <typename CRTP>
  struct MutableDAGViewBase
      : TargetView::StorageType::template MutableDAGViewBase<CRTP> {};

  MOVE_ONLY(OverlayDAGStorage);

  OverlayDAGStorage();
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
  struct NodeStorage : CombineBases<AllNodeFeatures> {};
  struct EdgeStorage : CombineBases<AllEdgeFeatures> {};

  auto GetTarget();
  auto GetTarget() const;

  std::decay_t<Target> target_ = {};
  std::unordered_map<NodeId, NodeStorage> node_storage_;
  std::unordered_map<NodeId, EdgeStorage> edge_storage_;
};
#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename CRTP, typename Tag>
template <typename F>
bool FeatureConstView<Overlay, CRTP, Tag>::IsOverlaid() const {
  auto& element_view = static_cast<const CRTP&>(*this);
  auto id = element_view.GetId();
  auto& storage = element_view.GetDAG().GetStorage();
  if constexpr (std::is_same_v<decltype(id), NodeId>) {
    if (storage.GetTarget().template ContainsId<CRTP>(id)) {
      auto it = tuple_get<OverlayFeatureStorageType<NodeId, F>, ContainerEquivalent>(
                    storage.replaced_node_storage_)
                    .find(id);
      return it != tuple_get<OverlayFeatureStorageType<NodeId, F>, ContainerEquivalent>(
                       storage.replaced_node_storage_)
                       .end();
    } else {
      return true;
    }
  } else {
    if (storage.GetTarget().template ContainsId<CRTP>(id)) {
      auto it = tuple_get<OverlayFeatureStorageType<EdgeId, F>, ContainerEquivalent>(
                    storage.replaced_edge_storage_)
                    .find(id);
      return it != tuple_get<OverlayFeatureStorageType<EdgeId, F>, ContainerEquivalent>(
                       storage.replaced_edge_storage_)
                       .end();
    } else {
      return true;
    }
  }
}

template <typename CRTP, typename Tag>
bool FeatureConstView<Overlay, CRTP, Tag>::IsAppended() const {
  auto& element_view = static_cast<const CRTP&>(*this);
  auto id = element_view.GetId();
  auto target_dag = element_view.GetDAG().GetStorage().GetTarget();
  return not target_dag.template ContainsId<CRTP>(id);
}

template <typename CRTP, typename Tag>
template <typename F>
auto FeatureMutableView<Overlay, CRTP, Tag>::SetOverlay() const {
  auto& element_view = static_cast<const CRTP&>(*this);
  auto id = element_view.GetId();
  const auto& storage = element_view.GetDAG().GetStorage().GetTargetStorage();

  static_assert(
      CRTP::template contains_feature<F>,
      "Attempted to SetOverlay on a Feature not supported by given DAG Element.");
  if constexpr (std::is_same_v<decltype(id), NodeId>) {
    Assert(storage.GetTarget().template ContainsId<CRTP>(id));
    auto& replaced_node_storage = tuple_get<OverlayFeatureStorageType<NodeId, F>,
                                            ContainerEquivalent>(
        element_view.GetDAG().GetStorage().GetTargetStorage().replaced_node_storage_);
    Assert(replaced_node_storage.find(id) == replaced_node_storage.end());
    if constexpr (std::is_copy_assignable_v<F>) {
      replaced_node_storage[id] =
          storage.GetTarget().template GetFeatureStorage<F>(id).get();
    } else {
      replaced_node_storage[id] =
          storage.GetTarget().template GetFeatureStorage<F>(id).get().Copy(
              static_cast<const CRTP*>(this));
    }
  } else {
    Assert(storage.GetTarget().template ContainsId<CRTP>(id));
    auto& replaced_edge_storage = tuple_get<OverlayFeatureStorageType<EdgeId, F>,
                                            ContainerEquivalent>(
        element_view.GetDAG().GetStorage().GetTargetStorage().replaced_edge_storage_);
    Assert(replaced_edge_storage.find(id) == replaced_edge_storage.end());
    if constexpr (std::is_copy_assignable_v<F>) {
      replaced_edge_storage[id] =
          storage.GetTarget().template GetFeatureStorage<F>(id).get();
    } else {
      replaced_edge_storage[id] =
          storage.GetTarget().template GetFeatureStorage<F>(id).get().Copy(
              static_cast<const CRTP*>(this));
    }
  }
  return element_view;
}

template <typename CRTP, typename Tag>
auto FeatureConstView<OverlayDAG, CRTP, Tag>::GetOriginal() const {
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetStorage().GetTargetStorage().GetTarget();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<OverlayDAG, CRTP, Tag>::HaveOverlays() const {
  auto& storage = static_cast<const CRTP&>(*this).GetStorage();

  return not(std::apply([](auto&&... ts) -> bool { return (ts.empty() and ...); },
                        storage.replaced_node_storage_) and
             std::apply([](auto&&... ts) -> bool { return (ts.empty() and ...); },
                        storage.replaced_edge_storage_) and
             storage.added_node_storage_.empty() and
             storage.added_edge_storage_.empty());
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <Component C, typename Feature>
inline constexpr bool
    OverlayDAGStorage<ShortName, Target, ViewBase>::contains_element_feature =
        TargetView::StorageType::template contains_element_feature<C, Feature>;

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <template <typename, typename> typename Base>
DAGView<typename OverlayDAGStorage<ShortName, Target, ViewBase>::Self, Base>
OverlayDAGStorage<ShortName, Target, ViewBase>::View() {
  return DAGView<Self, Base>{static_cast<Self&>(*this)};
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <template <typename, typename> typename Base>
DAGView<const typename OverlayDAGStorage<ShortName, Target, ViewBase>::Self, Base>
OverlayDAGStorage<ShortName, Target, ViewBase>::View() const {
  return DAGView<const Self, Base>{static_cast<const Self&>(*this)};
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
NodeId OverlayDAGStorage<ShortName, Target, ViewBase>::AppendNode() {
  auto result = GetNextAvailableId<Component::Node, TargetView>();
  added_node_storage_.push_back({});
  return result;
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
EdgeId OverlayDAGStorage<ShortName, Target, ViewBase>::AppendEdge() {
  auto result = GetNextAvailableId<Component::Edge, TargetView>();
  added_edge_storage_.push_back({});
  return result;
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
void OverlayDAGStorage<ShortName, Target, ViewBase>::AddNode(NodeId id) {
  View().Overlay(id);
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
void OverlayDAGStorage<ShortName, Target, ViewBase>::AddEdge(EdgeId id) {
  View().Overlay(id);
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <typename VT>
size_t OverlayDAGStorage<ShortName, Target, ViewBase>::GetNodesCount() const {
  return GetTarget().GetNodesCount() + added_node_storage_.size();
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <typename VT>
size_t OverlayDAGStorage<ShortName, Target, ViewBase>::GetEdgesCount() const {
  return GetTarget().GetEdgesCount() + added_edge_storage_.size();
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <typename VT>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetNodes() const {
  auto target_nodes = GetTarget().GetStorage().template GetNodes<VT>();
  auto first_added = GetTarget().template GetNextAvailableNodeId<VT>();
  auto added_nodes =
      ranges::views::indices(first_added.value,
                             first_added.value + added_node_storage_.size()) |
      Transform::ToId<Component::Node>();
  return ranges::views::concat(target_nodes, added_nodes);
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <typename VT>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetEdges() const {
  auto target_edges = GetTarget().GetStorage().template GetEdges<VT>();
  auto first_added = GetTarget().template GetNextAvailableEdgeId<VT>();
  auto added_edges =
      ranges::views::indices(first_added.value,
                             first_added.value + added_edge_storage_.size()) |
      Transform::ToId<Component::Edge>();
  return ranges::views::concat(target_edges, added_edges);
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
void OverlayDAGStorage<ShortName, Target, ViewBase>::InitializeNodes(size_t size) {
  if (size < GetTarget().GetNodesCount()) {
    Fail("Overlayed DAG can only be grown");
  }
  added_node_storage_.resize(size - GetTarget().GetNodesCount());
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
void OverlayDAGStorage<ShortName, Target, ViewBase>::InitializeEdges(size_t size) {
  if (size < GetTarget().GetEdgesCount()) {
    Fail("Overlayed DAG can only be grown");
  }
  added_edge_storage_.resize(size - GetTarget().GetEdgesCount());
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <typename F>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetFeatureStorage() {
  return GetTarget().template GetFeatureStorage<F>();
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <typename F>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetFeatureStorage() const {
  return GetTarget().template GetFeatureStorage<F>();
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <typename F>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetFeatureStorage(NodeId id) {
  return GetFeatureStorageImpl<F>(*this, id);
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <typename F>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetFeatureStorage(
    NodeId id) const {
  return GetFeatureStorageImpl<F>(*this, id);
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <typename F>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetFeatureStorage(EdgeId id) {
  return GetFeatureStorageImpl<F>(*this, id);
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <typename F>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetFeatureStorage(
    EdgeId id) const {
  return GetFeatureStorageImpl<F>(*this, id);
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <Component C, typename F>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetFeatureExtraStorage() {
  return GetTarget().template GetFeatureExtraStorage<C, F>();
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <Component C, typename F>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetFeatureExtraStorage() const {
  return GetTarget().template GetFeatureExtraStorage<C, F>();
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
OverlayDAGStorage<ShortName, Target, ViewBase>::OverlayDAGStorage(Target&& target)
    : target_{std::move(target)} {}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetTarget() {
  return ViewOf(target_);
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetTarget() const {
  return ViewOf(target_);
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <typename F, typename OverlayStorageType>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetFeatureStorageImpl(
    OverlayStorageType& self, NodeId id) {
  constexpr bool is_mutable = (not std::is_const_v<OverlayStorageType>) and
                              OverlayStorageType::TargetView::is_mutable;
  using Result =
      std::conditional_t<is_mutable, typename OverlayFeatureType<F>::mutable_view_type,
                         typename OverlayFeatureType<F>::const_view_type>;

  if (self.GetTarget().template ContainsId<TargetView>(id)) {
    auto it = tuple_get<OverlayFeatureStorageType<NodeId, F>, ContainerEquivalent>(
                  self.replaced_node_storage_)
                  .find(id);
    if (it == tuple_get<OverlayFeatureStorageType<NodeId, F>, ContainerEquivalent>(
                  self.replaced_node_storage_)
                  .end()) {
      if constexpr (is_mutable) {
        Fail("Can't modify non-overlaid node");
      } else {
        return Result{std::cref(self.GetTarget().template GetFeatureStorage<F>(id))};
      }
    } else {
      return Result{std::ref(it->second)};
    }
  } else {
    return Result{std::ref(tuple_get<typename OverlayFeatureType<F>::store_type,
                                     FeatureEquivalent>(self.added_node_storage_.at(
        id.value -
        self.GetTarget().template GetNextAvailableNodeId<TargetView>().value)))};
  }
}

template <typename ShortName, typename Target,
          template <typename, typename> typename ViewBase>
template <typename F, typename OverlayStorageType>
auto OverlayDAGStorage<ShortName, Target, ViewBase>::GetFeatureStorageImpl(
    OverlayStorageType& self, EdgeId id) {
  constexpr bool is_mutable = (not std::is_const_v<OverlayStorageType>) and
                              OverlayStorageType::TargetView::is_mutable;

  using Result =
      std::conditional_t<is_mutable, typename OverlayFeatureType<F>::mutable_view_type,
                         typename OverlayFeatureType<F>::const_view_type>;

  if (self.GetTarget().template ContainsId<TargetView>(id)) {
    auto it = tuple_get<OverlayFeatureStorageType<EdgeId, F>, ContainerEquivalent>(
                  self.replaced_edge_storage_)
                  .find(id);
    if (it == tuple_get<OverlayFeatureStorageType<EdgeId, F>, ContainerEquivalent>(
                  self.replaced_edge_storage_)
                  .end()) {
      if constexpr (is_mutable) {
        Fail("Can't modify non-overlaid edge");
      } else {
        const auto& result = self.GetTarget().template GetFeatureStorage<F>(id);
        return Result{result};
      }
    } else {
      return Result{std::ref(it->second)};
    }
  } else {
    return Result{std::ref(tuple_get<typename OverlayFeatureType<F>::store_type,
                                     FeatureEquivalent>(self.added_edge_storage_.at(
        id.value -
        self.GetTarget().template GetNextAvailableEdgeId<TargetView>().value)))};
  }
}

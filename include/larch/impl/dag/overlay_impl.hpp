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
    if (id.value < storage.GetTarget().GetNodesCount()) {
      auto it = std::get<std::unordered_map<NodeId, F>>(storage.replaced_node_storage_)
                    .find(id);
      return it !=
             std::get<std::unordered_map<NodeId, F>>(storage.replaced_node_storage_)
                 .end();
    } else {
      return true;
    }
  } else {
    if (id.value < storage.GetTarget().GetEdgesCount()) {
      auto it = std::get<std::unordered_map<EdgeId, F>>(storage.replaced_edge_storage_)
                    .find(id);
      return it !=
             std::get<std::unordered_map<EdgeId, F>>(storage.replaced_edge_storage_)
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
  if constexpr (std::is_same_v<decltype(id), NodeId>) {
    return id.value >= target_dag.GetNodesCount();
  } else {
    return id.value >= target_dag.GetEdgesCount();
  }
}

template <typename CRTP, typename Tag>
template <typename F>
auto FeatureMutableView<Overlay, CRTP, Tag>::SetOverlay() const {
  auto& element_view = static_cast<const CRTP&>(*this);
  auto id = element_view.GetId();
  auto& storage = element_view.GetDAG().GetStorage().GetTargetStorage();
  static_assert(
      CRTP::template contains_feature<F>,
      "Attempted to SetOverlay on a Feature not supported by given DAG Element.");
  if constexpr (std::is_same_v<decltype(id), NodeId>) {
    Assert(id.value < storage.GetTarget().GetNodesCount());
    auto& replaced_node_storage =
        std::get<std::unordered_map<NodeId, F>>(storage.replaced_node_storage_);
    Assert(replaced_node_storage.find(id) == replaced_node_storage.end());
    if constexpr (std::is_copy_assignable_v<F>) {
      replaced_node_storage[id] = storage.GetTarget().template GetFeatureStorage<F>(id);
    } else {
      replaced_node_storage[id] =
          storage.GetTarget().template GetFeatureStorage<F>(id).Copy();
    }
  } else {
    Assert(id.value < storage.GetTarget().GetEdgesCount());
    auto& replaced_edge_storage =
        std::get<std::unordered_map<EdgeId, F>>(storage.replaced_edge_storage_);
    Assert(replaced_edge_storage.find(id) == replaced_edge_storage.end());
    if constexpr (std::is_copy_assignable_v<F>) {
      replaced_edge_storage[id] = storage.GetTarget().template GetFeatureStorage<F>(id);
    } else {
      replaced_edge_storage[id] =
          storage.GetTarget().template GetFeatureStorage<F>(id).Copy();
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
  return not(storage.replaced_node_storage_.empty() and
             storage.replaced_edge_storage_.empty() and
             storage.added_node_storage_.empty() and
             storage.added_edge_storage_.empty());
}

template <typename ShortName, typename Target>
template <Component C, typename Feature>
inline constexpr bool OverlayDAGStorage<ShortName, Target>::contains_element_feature =
    TargetView::StorageType::template contains_element_feature<C, Feature>;

template <typename ShortName, typename Target>
auto OverlayDAGStorage<ShortName, Target>::View() {
  return DAGView<Self>{static_cast<Self&>(*this)};
}

template <typename ShortName, typename Target>
auto OverlayDAGStorage<ShortName, Target>::View() const {
  return DAGView<const Self>{static_cast<const Self&>(*this)};
}

template <typename ShortName, typename Target>
NodeId OverlayDAGStorage<ShortName, Target>::AppendNode() {
  added_node_storage_.push_back({});
  return {GetNodesCount() - 1};
}

template <typename ShortName, typename Target>
EdgeId OverlayDAGStorage<ShortName, Target>::AppendEdge() {
  added_edge_storage_.push_back({});
  return {GetEdgesCount() - 1};
}

template <typename ShortName, typename Target>
void OverlayDAGStorage<ShortName, Target>::AddNode(NodeId id) {
  View().Overlay(id);
}

template <typename ShortName, typename Target>
void OverlayDAGStorage<ShortName, Target>::AddEdge(EdgeId id) {
  View().Overlay(id);
}

template <typename ShortName, typename Target>
size_t OverlayDAGStorage<ShortName, Target>::GetNodesCount() const {
  return GetTarget().GetNodesCount() + added_node_storage_.size();
}

template <typename ShortName, typename Target>
size_t OverlayDAGStorage<ShortName, Target>::GetEdgesCount() const {
  return GetTarget().GetEdgesCount() + added_edge_storage_.size();
}

template <typename ShortName, typename Target>
auto OverlayDAGStorage<ShortName, Target>::GetNodes() const {
  return ranges::views::indices(GetNodesCount()) |
         ranges::views::transform([this](size_t i) {
           return ElementView{this->View(), NodeId{i}};
         });
}

template <typename ShortName, typename Target>
auto OverlayDAGStorage<ShortName, Target>::GetEdges() const {
  return ranges::views::indices(GetEdgesCount()) |
         ranges::views::transform([this](size_t i) {
           return ElementView{this->View(), EdgeId{i}};
         });
}

template <typename ShortName, typename Target>
void OverlayDAGStorage<ShortName, Target>::InitializeNodes(size_t size) {
  if (size < GetTarget().GetNodesCount()) {
    Fail("Overlayed DAG can only be grown");
  }
  added_node_storage_.resize(size - GetTarget().GetNodesCount());
}

template <typename ShortName, typename Target>
void OverlayDAGStorage<ShortName, Target>::InitializeEdges(size_t size) {
  if (size < GetTarget().GetEdgesCount()) {
    Fail("Overlayed DAG can only be grown");
  }
  added_edge_storage_.resize(size - GetTarget().GetEdgesCount());
}

template <typename ShortName, typename Target>
template <typename F>
auto& OverlayDAGStorage<ShortName, Target>::GetFeatureStorage() {
  return GetTarget().template GetFeatureStorage<F>();
}

template <typename ShortName, typename Target>
template <typename F>
const auto& OverlayDAGStorage<ShortName, Target>::GetFeatureStorage() const {
  return GetTarget().template GetFeatureStorage<F>();
}

template <typename ShortName, typename Target>
template <typename F>
auto& OverlayDAGStorage<ShortName, Target>::GetFeatureStorage(NodeId id) {
  return GetFeatureStorageImpl<F>(*this, id);
}

template <typename ShortName, typename Target>
template <typename F>
const auto& OverlayDAGStorage<ShortName, Target>::GetFeatureStorage(NodeId id) const {
  return GetFeatureStorageImpl<F>(*this, id);
}

template <typename ShortName, typename Target>
template <typename F>
auto& OverlayDAGStorage<ShortName, Target>::GetFeatureStorage(EdgeId id) {
  return GetFeatureStorageImpl<F>(*this, id);
}

template <typename ShortName, typename Target>
template <typename F>
const auto& OverlayDAGStorage<ShortName, Target>::GetFeatureStorage(EdgeId id) const {
  return GetFeatureStorageImpl<F>(*this, id);
}

template <typename ShortName, typename Target>
template <Component C, typename F>
auto& OverlayDAGStorage<ShortName, Target>::GetFeatureExtraStorage() {
  return GetTarget().template GetFeatureExtraStorage<C, F>();
}

template <typename ShortName, typename Target>
template <Component C, typename F>
const auto& OverlayDAGStorage<ShortName, Target>::GetFeatureExtraStorage() const {
  return GetTarget().template GetFeatureExtraStorage<C, F>();
}

template <typename ShortName, typename Target>
OverlayDAGStorage<ShortName, Target>::OverlayDAGStorage(Target&& target)
    : target_{std::move(target)} {}

template <typename ShortName, typename Target>
auto OverlayDAGStorage<ShortName, Target>::GetTarget() {
  return ViewOf(target_);
}

template <typename ShortName, typename Target>
auto OverlayDAGStorage<ShortName, Target>::GetTarget() const {
  return ViewOf(target_);
}

template <typename ShortName, typename Target>
template <typename F, typename OverlayStorageType>
auto OverlayDAGStorage<ShortName, Target>::GetFeatureStorageImpl(
    OverlayStorageType& self, NodeId id)
    -> std::conditional_t<not std::is_const_v<OverlayStorageType> and
                              OverlayStorageType::TargetView::is_mutable,
                          F&, const F&> {
  Assert(id.value < self.GetNodesCount());
  if (id.value < self.GetTarget().GetNodesCount()) {
    auto it =
        std::get<std::unordered_map<NodeId, F>>(self.replaced_node_storage_).find(id);
    if (it ==
        std::get<std::unordered_map<NodeId, F>>(self.replaced_node_storage_).end()) {
      if constexpr (not std::is_const_v<OverlayStorageType> and
                    OverlayStorageType::TargetView::is_mutable) {
        Fail("Can't modify non-overlaid node");
      }
      return self.GetTarget().template GetFeatureStorage<F>(id);
    } else {
      return it->second;
    }
  } else {
    return std::get<F>(
        self.added_node_storage_.at(id.value - self.GetTarget().GetNodesCount()));
  }
}

template <typename ShortName, typename Target>
template <typename F, typename OverlayStorageType>
auto OverlayDAGStorage<ShortName, Target>::GetFeatureStorageImpl(
    OverlayStorageType& self, EdgeId id)
    -> std::conditional_t<not std::is_const_v<OverlayStorageType> and
                              OverlayStorageType::TargetView::is_mutable,
                          F&, const F&> {
  if (id.value < self.GetTarget().GetEdgesCount()) {
    auto it =
        std::get<std::unordered_map<EdgeId, F>>(self.replaced_edge_storage_).find(id);
    if (it ==
        std::get<std::unordered_map<EdgeId, F>>(self.replaced_edge_storage_).end()) {
      if constexpr (not std::is_const_v<OverlayStorageType> and
                    OverlayStorageType::TargetView::is_mutable) {
        Fail("Can't modify non-overlaid edge");
      }
      return self.GetTarget().template GetFeatureStorage<F>(id);
    } else {
      return it->second;
    }
  } else {
    return std::get<F>(
        self.added_edge_storage_.at(id.value - self.GetTarget().GetEdgesCount()));
  }
}

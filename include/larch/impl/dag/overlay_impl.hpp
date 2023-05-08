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
  auto& storage = element_view.GetDAG().GetStorage();
  static_assert(
      CRTP::template contains_feature<F>,
      "Attempted to SetOverlay on a Feature not supported by given DAG Element.");
  if constexpr (std::is_same_v<decltype(id), NodeId>) {
    if (id.value < storage.GetTarget().GetNodesCount()) {
      if constexpr (std::is_copy_assignable_v<F>) {
        std::get<std::unordered_map<NodeId, F>>(storage.replaced_node_storage_)[id] =
            storage.GetTarget().template GetFeatureStorage<F>(id);
      } else {
        std::get<std::unordered_map<NodeId, F>>(storage.replaced_node_storage_)[id] =
            storage.GetTarget().template GetFeatureStorage<F>(id).Copy();
      }
    } else {
      std::ignore = GetOrInsert(storage.added_node_storage_, id);
    }
  } else {
    if (id.value < storage.GetTarget().GetEdgesCount()) {
      if constexpr (std::is_copy_assignable_v<F>) {
        std::get<std::unordered_map<EdgeId, F>>(storage.replaced_edge_storage_)[id] =
            storage.GetTarget().template GetFeatureStorage<F>(id);
      } else {
        std::get<std::unordered_map<EdgeId, F>>(storage.replaced_edge_storage_)[id] =
            storage.GetTarget().template GetFeatureStorage<F>(id).Copy();
      }
    } else {
      std::ignore = GetOrInsert(storage.added_edge_storage_, id);
    }
  }
  return element_view;
}

template <typename CRTP, typename Tag>
auto FeatureConstView<OverlayDAG, CRTP, Tag>::GetOriginal() const {
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetStorage().GetTarget();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<OverlayDAG, CRTP, Tag>::HaveOverlays() const {
  auto& storage = static_cast<const CRTP&>(*this).GetStorage();
  return not(storage.replaced_node_storage_.empty() and
             storage.replaced_edge_storage_.empty() and
             storage.added_node_storage_.empty() and
             storage.added_edge_storage_.empty());
}

template <typename Target>
template <typename Id, typename Feature>
inline constexpr bool OverlayDAGStorage<Target>::contains_element_feature =
    TargetView::StorageType::template contains_element_feature<Id, Feature>;

template <typename Target>
OverlayDAGStorage<Target>::OverlayDAGStorage(Target&& target)
    : target_{std::forward<Target>(target)} {}

template <typename Target>
auto OverlayDAGStorage<Target>::View() {
  return DAGView<OverlayDAGStorage<Target>>{*this};
}

template <typename Target>
auto OverlayDAGStorage<Target>::View() const {
  return DAGView<const OverlayDAGStorage<Target>>{*this};
}

template <typename Target>
NodeId OverlayDAGStorage<Target>::AppendNode() {
  added_node_storage_.push_back({});
  return {GetNodesCount() - 1};
}

template <typename Target>
EdgeId OverlayDAGStorage<Target>::AppendEdge() {
  added_edge_storage_.push_back({});
  return {GetEdgesCount() - 1};
}

template <typename Target>
void OverlayDAGStorage<Target>::AddNode(NodeId id) {
  View().Overlay(id);
}

template <typename Target>
void OverlayDAGStorage<Target>::AddEdge(EdgeId id) {
  View().Overlay(id);
}

template <typename Target>
size_t OverlayDAGStorage<Target>::GetNodesCount() const {
  return GetTarget().GetNodesCount() + added_node_storage_.size();
}

template <typename Target>
size_t OverlayDAGStorage<Target>::GetEdgesCount() const {
  return GetTarget().GetEdgesCount() + added_edge_storage_.size();
}

template <typename Target>
auto OverlayDAGStorage<Target>::GetNodes() const {
  return ranges::views::indices(GetNodesCount()) |
         ranges::views::transform([this](size_t i) {
           return ElementView{this->View(), NodeId{i}};
         });
}

template <typename Target>
auto OverlayDAGStorage<Target>::GetEdges() const {
  return ranges::views::indices(GetEdgesCount()) |
         ranges::views::transform([this](size_t i) {
           return ElementView{this->View(), EdgeId{i}};
         });
}

template <typename Target>
void OverlayDAGStorage<Target>::InitializeNodes(size_t size) {
  if (size < GetTarget().GetNodesCount()) {
    Fail("Overlayed DAG can only be grown");
  }
  added_node_storage_.resize(size - GetTarget().GetNodesCount());
}

template <typename Target>
template <typename F>
auto& OverlayDAGStorage<Target>::GetFeatureStorage() {
  return GetTarget().template GetFeatureStorage<F>();
}

template <typename Target>
template <typename F>
const auto& OverlayDAGStorage<Target>::GetFeatureStorage() const {
  return GetTarget().template GetFeatureStorage<F>();
}

template <typename Target>
template <typename F>
auto& OverlayDAGStorage<Target>::GetFeatureStorage(NodeId id) {
  return GetFeatureStorageImpl<F>(*this, id);
}

template <typename Target>
template <typename F>
const auto& OverlayDAGStorage<Target>::GetFeatureStorage(NodeId id) const {
  return GetFeatureStorageImpl<F>(*this, id);
}

template <typename Target>
template <typename F>
auto& OverlayDAGStorage<Target>::GetFeatureStorage(EdgeId id) {
  return GetFeatureStorageImpl<F>(*this, id);
}

template <typename Target>
template <typename F>
const auto& OverlayDAGStorage<Target>::GetFeatureStorage(EdgeId id) const {
  return GetFeatureStorageImpl<F>(*this, id);
}

template <typename Target>
template <typename Id, typename F>
auto& OverlayDAGStorage<Target>::GetFeatureExtraStorage() {
  return GetTarget().template GetFeatureExtraStorage<Id, F>();
}

template <typename Target>
template <typename Id, typename F>
const auto& OverlayDAGStorage<Target>::GetFeatureExtraStorage() const {
  return GetTarget().template GetFeatureExtraStorage<Id, F>();
}

template <typename Target>
auto OverlayDAGStorage<Target>::GetTarget() {
  return ViewOf(target_);
}

template <typename Target>
auto OverlayDAGStorage<Target>::GetTarget() const {
  return ViewOf(target_);
}

template <typename Target>
template <typename F, typename OverlayStorageType>
auto OverlayDAGStorage<Target>::GetFeatureStorageImpl(OverlayStorageType& self,
                                                      NodeId id)
    -> std::conditional_t<not std::is_const_v<OverlayStorageType> and
                              OverlayStorageType::TargetView::is_mutable,
                          F&, const F&> {
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

template <typename Target>
template <typename F, typename OverlayStorageType>
auto OverlayDAGStorage<Target>::GetFeatureStorageImpl(OverlayStorageType& self,
                                                      EdgeId id)
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

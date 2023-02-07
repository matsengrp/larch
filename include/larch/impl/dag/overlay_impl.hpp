#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename CRTP, typename Tag>
bool FeatureConstView<Overlay, CRTP, Tag>::IsOverlaid() const {
  auto& element_view = static_cast<const CRTP&>(*this);
  auto id = element_view.GetId();
  auto& storage = element_view.GetDAG().GetStorage();
  if constexpr (std::is_same_v<decltype(id), NodeId>) {
    if (id.value < storage.GetTarget().GetNodesCount()) {
      auto it = storage.replaced_node_storage_.find(id);
      return (it != storage.replaced_node_storage_.end());
    } else {
      return true;
    }
  } else {
    if (id.value < storage.GetTarget().GetEdgesCount()) {
      auto it = storage.replaced_edge_storage_.find(id);
      return (it != storage.replaced_edge_storage_.end());
    } else {
      return true;
    }
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Overlay, CRTP, Tag>::Overlay() {
  auto& element_view = static_cast<const CRTP&>(*this);
  auto id = element_view.GetId();
  auto& storage = element_view.GetDAG().GetStorage();
  if constexpr (std::is_same_v<decltype(id), NodeId>) {
    if (id.value < storage.GetTarget().GetNodesCount()) {
      storage.replaced_node_storage_[id] = {};
    } else {
      std::ignore = GetOrInsert(storage.added_node_storage_, id);
    }
  } else {
    if (id.value < storage.GetTarget().GetEdgesCount()) {
      storage.replaced_edge_storage_[id] = {};
    } else {
      std::ignore = GetOrInsert(storage.replaced_edge_storage_, id);
    }
  }
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
  return {GetNodesCount() - 1};
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
    throw std::runtime_error("Overlayed DAG can only be grown");
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
  if (id.value < GetTarget().GetNodesCount()) {
    auto it = replaced_node_storage_.find(id);
    if (it == replaced_node_storage_.end()) {
      return GetTarget().template GetFeatureStorage<F>(id);
    } else {
      return std::get<F>(it->second);
    }
  } else {
    return std::get<F>(added_node_storage_.at(id.value - GetTarget().GetNodesCount()));
  }
}

template <typename Target>
template <typename F>
const auto& OverlayDAGStorage<Target>::GetFeatureStorage(NodeId id) const {
  if (id.value < GetTarget().GetNodesCount()) {
    auto it = replaced_node_storage_.find(id);
    if (it == replaced_node_storage_.end()) {
      return GetTarget().template GetFeatureStorage<F>(id);
    } else {
      return std::get<F>(it->second);
    }
  } else {
    return std::get<F>(added_node_storage_.at(id.value - GetTarget().GetNodesCount()));
  }
}

template <typename Target>
template <typename F>
auto& OverlayDAGStorage<Target>::GetFeatureStorage(EdgeId id) {
  if (id.value < GetTarget().GetEdgesCount()) {
    auto it = replaced_edge_storage_.find(id);
    if (it == replaced_edge_storage_.end()) {
      return GetTarget().template GetFeatureStorage<F>(id);
    } else {
      return std::get<F>(it->second);
    }
  } else {
    return std::get<F>(added_edge_storage_.at(id.value - GetTarget().GetEdgesCount()));
  }
}

template <typename Target>
template <typename F>
const auto& OverlayDAGStorage<Target>::GetFeatureStorage(EdgeId id) const {
  if (id.value < GetTarget().GetEdgesCount()) {
    auto it = replaced_edge_storage_.find(id);
    if (it == replaced_edge_storage_.end()) {
      return GetTarget().template GetFeatureStorage<F>(id);
    } else {
      return std::get<F>(it->second);
    }
  } else {
    return std::get<F>(added_edge_storage_.at(id.value - GetTarget().GetEdgesCount()));
  }
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
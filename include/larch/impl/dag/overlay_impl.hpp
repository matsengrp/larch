#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Target>
template <typename Id, typename Feature>
inline constexpr bool OverlayDAGStorage<Target>::contains_element_feature =
    TargetView::StorageType::template contains_element_feature<Id, Feature>;

template <typename Target>
OverlayDAGStorage<Target>::OverlayDAGStorage() = default;

template <typename Target>
OverlayDAGStorage<Target>::OverlayDAGStorage(Target&& target)
    : target_{std::move(target)} {}

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
  return GetTarget().AppendNode().GetId();
}

template <typename Target>
EdgeId OverlayDAGStorage<Target>::AppendEdge() {
  return GetTarget().AppendEdge().GetId();
}

template <typename Target>
void OverlayDAGStorage<Target>::AddNode(NodeId id) {
  GetTarget().AddNode(id);
}

template <typename Target>
void OverlayDAGStorage<Target>::AddEdge(EdgeId id) {
  GetTarget().AddEdge(id);
}

template <typename Target>
size_t OverlayDAGStorage<Target>::GetNodesCount() const {
  return GetTarget().GetNodesCount();
}

template <typename Target>
size_t OverlayDAGStorage<Target>::GetEdgesCount() const {
  return GetTarget().GetEdgesCount();
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
  GetTarget().InitializeNodes(size);
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
  auto it = node_storage_.find(id);
  if (it == node_storage_.end()) {
    return GetTarget().template GetFeatureStorage<F>(id);
  } else {
    return std::get<F>(it->second);
  }
}

template <typename Target>
template <typename F>
const auto& OverlayDAGStorage<Target>::GetFeatureStorage(NodeId id) const {
  auto it = node_storage_.find(id);
  if (it == node_storage_.end()) {
    return GetTarget().template GetFeatureStorage<F>(id);
  } else {
    return std::get<F>(it->second);
  }
}

template <typename Target>
template <typename F>
auto& OverlayDAGStorage<Target>::GetFeatureStorage(EdgeId id) {
  auto it = edge_storage_.find(id);
  if (it == edge_storage_.end()) {
    return GetTarget().template GetFeatureStorage<F>(id);
  } else {
    return std::get<F>(it->second);
  }
}

template <typename Target>
template <typename F>
const auto& OverlayDAGStorage<Target>::GetFeatureStorage(EdgeId id) const {
  auto it = edge_storage_.find(id);
  if (it == edge_storage_.end()) {
    return GetTarget().template GetFeatureStorage<F>(id);
  } else {
    return std::get<F>(it->second);
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
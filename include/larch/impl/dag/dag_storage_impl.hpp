#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename NC, typename EC, typename... Fs>
template <typename Id, typename Feature>
inline constexpr bool DAGStorage<NC, EC, Fs...>::contains_element_feature = [] {
  if constexpr (std::is_same_v<Id, NodeId>) {
    return NC::template contains_element_feature<Feature>;
  } else {
    return EC::template contains_element_feature<Feature>;
  }
}();

template <typename NC, typename EC, typename... Fs>
auto DAGStorage<NC, EC, Fs...>::View() {
  return DAGView<DAGStorage<NC, EC, Fs...>>{*this};
}

template <typename NC, typename EC, typename... Fs>
auto DAGStorage<NC, EC, Fs...>::View() const {
  return DAGView<const DAGStorage<NC, EC, Fs...>>{*this};
}

template <typename NC, typename EC, typename... Fs>
NodeId DAGStorage<NC, EC, Fs...>::AppendNode() {
  return nodes_container_.Append();
}

template <typename NC, typename EC, typename... Fs>
EdgeId DAGStorage<NC, EC, Fs...>::AppendEdge() {
  return edges_container_.Append();
}

template <typename NC, typename EC, typename... Fs>
void DAGStorage<NC, EC, Fs...>::AddNode(NodeId id) {
  nodes_container_.Add(id);
}

template <typename NC, typename EC, typename... Fs>
void DAGStorage<NC, EC, Fs...>::AddEdge(EdgeId id) {
  edges_container_.Add(id);
}

template <typename NC, typename EC, typename... Fs>
size_t DAGStorage<NC, EC, Fs...>::GetNodesCount() const {
  return nodes_container_.GetCount();
}

template <typename NC, typename EC, typename... Fs>
size_t DAGStorage<NC, EC, Fs...>::GetEdgesCount() const {
  return edges_container_.GetCount();
}

template <typename NC, typename EC, typename... Fs>
void DAGStorage<NC, EC, Fs...>::InitializeNodes(size_t size) {
  nodes_container_.Initialize(size);
}

template <typename NC, typename EC, typename... Fs>
template <typename F>
auto& DAGStorage<NC, EC, Fs...>::GetFeatureStorage(NodeId id) {
  return nodes_container_.template GetFeatureStorage<F>(id);
}

template <typename NC, typename EC, typename... Fs>
template <typename F>
const auto& DAGStorage<NC, EC, Fs...>::GetFeatureStorage(NodeId id) const {
  return nodes_container_.template GetFeatureStorage<F>(id);
}

template <typename NC, typename EC, typename... Fs>
template <typename F>
auto& DAGStorage<NC, EC, Fs...>::GetFeatureStorage(EdgeId id) {
  return edges_container_.template GetFeatureStorage<F>(id);
}

template <typename NC, typename EC, typename... Fs>
template <typename F>
const auto& DAGStorage<NC, EC, Fs...>::GetFeatureStorage(EdgeId id) const {
  return edges_container_.template GetFeatureStorage<F>(id);
}

template <typename NC, typename EC, typename... Fs>
template <typename Id, typename F>
auto& DAGStorage<NC, EC, Fs...>::GetFeatureExtraStorage() {
  if constexpr (std::is_same_v<Id, NodeId>) {
    return nodes_container_.template GetFeatureExtraStorage<F>();
  } else {
    return edges_container_.template GetFeatureExtraStorage<F>();
  }
}

template <typename NC, typename EC, typename... Fs>
template <typename Id, typename F>
const auto& DAGStorage<NC, EC, Fs...>::GetFeatureExtraStorage() const {
  if constexpr (std::is_same_v<Id, NodeId>) {
    return nodes_container_.template GetFeatureExtraStorage<F>();
  } else {
    return edges_container_.template GetFeatureExtraStorage<F>();
  }
}

template <typename NC, typename EC, typename... Fs>
template <typename F>
auto& DAGStorage<NC, EC, Fs...>::GetFeatureStorage() {
  return std::get<F>(features_storage_);
}

template <typename NC, typename EC, typename... Fs>
template <typename F>
const auto& DAGStorage<NC, EC, Fs...>::GetFeatureStorage() const {
  return std::get<F>(features_storage_);
}

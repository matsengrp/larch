#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
template <typename Id, typename Feature>
inline constexpr bool DAGStorage<NodesContainerT, EdgesContainerT,
                                 Features...>::contains_element_feature = [] {
  if constexpr (std::is_same_v<Id, NodeId>) {
    return NodesContainerT::template contains_element_feature<Feature>;
  } else {
    return EdgesContainerT::template contains_element_feature<Feature>;
  }
}();

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
auto DAGStorage<NodesContainerT, EdgesContainerT, Features...>::View() {
  return DAGView<DAGStorage<NodesContainerT, EdgesContainerT, Features...>>{*this};
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
auto DAGStorage<NodesContainerT, EdgesContainerT, Features...>::View() const {
  return DAGView<const DAGStorage<NodesContainerT, EdgesContainerT, Features...>>{
      *this};
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
NodeId DAGStorage<NodesContainerT, EdgesContainerT, Features...>::AppendNode() {
  return nodes_container_.Append();
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
EdgeId DAGStorage<NodesContainerT, EdgesContainerT, Features...>::AppendEdge() {
  return edges_container_.Append();
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
void DAGStorage<NodesContainerT, EdgesContainerT, Features...>::AddNode(NodeId id) {
  nodes_container_.Add(id);
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
void DAGStorage<NodesContainerT, EdgesContainerT, Features...>::AddEdge(EdgeId id) {
  edges_container_.Add(id);
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
size_t DAGStorage<NodesContainerT, EdgesContainerT, Features...>::GetNodesCount()
    const {
  return nodes_container_.GetCount();
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
size_t DAGStorage<NodesContainerT, EdgesContainerT, Features...>::GetEdgesCount()
    const {
  return edges_container_.GetCount();
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
const auto& DAGStorage<NodesContainerT, EdgesContainerT, Features...>::GetNodes()
    const {
  return nodes_container_.elements_storage_;
}  // TODO

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
const auto& DAGStorage<NodesContainerT, EdgesContainerT, Features...>::GetEdges()
    const {
  return edges_container_.elements_storage_;
}  // TODO

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
void DAGStorage<NodesContainerT, EdgesContainerT, Features...>::InitializeNodes(
    size_t size) {
  nodes_container_.Initialize(size);
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
template <typename Feature>
auto& DAGStorage<NodesContainerT, EdgesContainerT, Features...>::GetFeatureStorage(
    NodeId id) {
  return nodes_container_.template GetFeatureStorage<Feature>(id);
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
template <typename Feature>
const auto& DAGStorage<NodesContainerT, EdgesContainerT,
                       Features...>::GetFeatureStorage(NodeId id) const {
  return nodes_container_.template GetFeatureStorage<Feature>(id);
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
template <typename Feature>
auto& DAGStorage<NodesContainerT, EdgesContainerT, Features...>::GetFeatureStorage(
    EdgeId id) {
  return edges_container_.template GetFeatureStorage<Feature>(id);
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
template <typename Feature>
const auto& DAGStorage<NodesContainerT, EdgesContainerT,
                       Features...>::GetFeatureStorage(EdgeId id) const {
  return edges_container_.template GetFeatureStorage<Feature>(id);
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
template <typename Id, typename Feature>
auto& DAGStorage<NodesContainerT, EdgesContainerT,
                 Features...>::GetFeatureExtraStorage() {
  if constexpr (std::is_same_v<Id, NodeId>) {
    return nodes_container_.template GetFeatureExtraStorage<Feature>();
  } else {
    return edges_container_.template GetFeatureExtraStorage<Feature>();
  }
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
template <typename Id, typename Feature>
const auto& DAGStorage<NodesContainerT, EdgesContainerT,
                       Features...>::GetFeatureExtraStorage() const {
  if constexpr (std::is_same_v<Id, NodeId>) {
    return nodes_container_.template GetFeatureExtraStorage<Feature>();
  } else {
    return edges_container_.template GetFeatureExtraStorage<Feature>();
  }
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
template <typename Feature>
auto& DAGStorage<NodesContainerT, EdgesContainerT, Features...>::GetFeatureStorage() {
  return std::get<Feature>(features_storage_);
}

template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
template <typename Feature>
const auto&
DAGStorage<NodesContainerT, EdgesContainerT, Features...>::GetFeatureStorage() const {
  return std::get<Feature>(features_storage_);
}

#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::DAGStorage(
    NodesContainerT&& nodes_container, EdgesContainerT&& edges_container,
    ExtraStorageT&& features_storage)
    : nodes_container_{std::move(nodes_container)},
      edges_container_{std::move(edges_container)},
      features_storage_{std::move(features_storage)} {}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
template <Component C, typename Feature>
inline constexpr bool DAGStorage<NodesContainerT, EdgesContainerT,
                                 ExtraStorageT>::contains_element_feature = [] {
  // NOLINTBEGIN
  if constexpr (C == Component::Node) {
    return NodesContainerT::template contains_element_feature<Feature>;
  } else {
    return EdgesContainerT::template contains_element_feature<Feature>;
  }
  // NOLINTEND
}();

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
DAGView<DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>>
DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::View() {
  return ViewType{*this};
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
DAGView<const DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>>
DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::View() const {
  return ConstViewType{*this};
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
NodeId DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::AppendNode() {
  return nodes_container_.Append();
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
EdgeId DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::AppendEdge() {
  return edges_container_.Append();
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
void DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::AddNode(NodeId id) {
  nodes_container_.Add(id);
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
void DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::AddEdge(EdgeId id) {
  edges_container_.Add(id);
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
size_t DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::GetNodesCount()
    const {
  return nodes_container_.GetCount();
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
size_t DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::GetEdgesCount()
    const {
  return edges_container_.GetCount();
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
void DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::InitializeNodes(
    size_t size) {
  nodes_container_.Initialize(size);
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
void DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::InitializeEdges(
    size_t size) {
  edges_container_.Initialize(size);
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
void DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::ClearNodes() {
  nodes_container_.Clear();
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
void DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::ClearEdges() {
  edges_container_.Clear();
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
template <typename Feature>
auto& DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::GetFeatureStorage(
    NodeId id) {
  return nodes_container_.template GetFeatureStorage<Feature>(id);
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
template <typename Feature>
const auto& DAGStorage<NodesContainerT, EdgesContainerT,
                       ExtraStorageT>::GetFeatureStorage(NodeId id) const {
  return nodes_container_.template GetFeatureStorage<Feature>(id);
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
template <typename Feature>
auto& DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::GetFeatureStorage(
    EdgeId id) {
  return edges_container_.template GetFeatureStorage<Feature>(id);
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
template <typename Feature>
const auto& DAGStorage<NodesContainerT, EdgesContainerT,
                       ExtraStorageT>::GetFeatureStorage(EdgeId id) const {
  return edges_container_.template GetFeatureStorage<Feature>(id);
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
template <Component C, typename Feature>
auto& DAGStorage<NodesContainerT, EdgesContainerT,
                 ExtraStorageT>::GetFeatureExtraStorage() {
  if constexpr (C == Component::Node) {
    return nodes_container_.template GetFeatureExtraStorage<Feature>();
  } else {
    return edges_container_.template GetFeatureExtraStorage<Feature>();
  }
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
template <Component C, typename Feature>
const auto& DAGStorage<NodesContainerT, EdgesContainerT,
                       ExtraStorageT>::GetFeatureExtraStorage() const {
  if constexpr (C == Component::Node) {
    return nodes_container_.template GetFeatureExtraStorage<Feature>();
  } else {
    return edges_container_.template GetFeatureExtraStorage<Feature>();
  }
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
template <typename Feature>
auto& DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::GetFeatureStorage() {
  return features_storage_.template GetFeatureStorage<Feature>();
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
template <typename Feature>
const auto&
DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::GetFeatureStorage() const {
  return features_storage_.template GetFeatureStorage<Feature>();
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
template <Component C>
auto& DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::GetContainer() {
  if constexpr (C == Component::Node) {
    return nodes_container_;
  } else {
    return edges_container_;
  }
}

template <typename NodesContainerT, typename EdgesContainerT, typename ExtraStorageT>
template <Component C>
const auto& DAGStorage<NodesContainerT, EdgesContainerT, ExtraStorageT>::GetContainer()
    const {
  if constexpr (C == Component::Node) {
    return nodes_container_;
  } else {
    return edges_container_;
  }
}

struct DefaultDAGStorage
    : DAGStorage<ElementsContainer<Component::Node, ElementStorage<Neighbors>>,
                 ElementsContainer<Component::Edge, ElementStorage<Endpoints>>,
                 ExtraStorage<Connections>> {};

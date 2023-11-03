#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT>::DAGStorage(
    NodesContainerT&& nodes_container, EdgesContainerT&& edges_container,
    ExtraStorageT&& features_storage)
    : nodes_container_{std::move(nodes_container)},
      edges_container_{std::move(edges_container)},
      features_storage_{std::move(features_storage)} {}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
template <Component C, typename Feature>
inline constexpr bool DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                                 ExtraStorageT>::contains_element_feature = [] {
  // NOLINTBEGIN
  if constexpr (C == Component::Node) {
    return NodesContainerT::template contains_element_feature<Feature>;
  } else {
    return EdgesContainerT::template contains_element_feature<Feature>;
  }
  // NOLINTEND
}();

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
typename DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                    ExtraStorageT>::ViewType
DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT>::View() {
  return DAGView<Self>{static_cast<Self&>(*this)};
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
typename DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                    ExtraStorageT>::ConstViewType
DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT>::View() const {
  return DAGView<const Self>{static_cast<const Self&>(*this)};
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
NodeId
DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT>::AppendNode() {
  return nodes_container_.Append();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
EdgeId
DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT>::AppendEdge() {
  return edges_container_.Append();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
void DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT>::AddNode(
    NodeId id) {
  nodes_container_.Add(id);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
void DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT>::AddEdge(
    EdgeId id) {
  edges_container_.Add(id);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
size_t DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                  ExtraStorageT>::GetNodesCount() const {
  return nodes_container_.GetCount();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
size_t DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                  ExtraStorageT>::GetEdgesCount() const {
  return edges_container_.GetCount();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
void DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                ExtraStorageT>::InitializeNodes(size_t size) {
  nodes_container_.Initialize(size);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
void DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                ExtraStorageT>::InitializeEdges(size_t size) {
  edges_container_.Initialize(size);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
void DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                ExtraStorageT>::ClearNodes() {
  nodes_container_.Clear();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
void DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                ExtraStorageT>::ClearEdges() {
  edges_container_.Clear();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
template <typename Feature>
auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                 ExtraStorageT>::GetFeatureStorage(NodeId id) {
  return nodes_container_.template GetFeatureStorage<Feature>(id);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
template <typename Feature>
const auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                       ExtraStorageT>::GetFeatureStorage(NodeId id) const {
  return nodes_container_.template GetFeatureStorage<Feature>(id);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
template <typename Feature>
auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                 ExtraStorageT>::GetFeatureStorage(EdgeId id) {
  return edges_container_.template GetFeatureStorage<Feature>(id);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
template <typename Feature>
const auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                       ExtraStorageT>::GetFeatureStorage(EdgeId id) const {
  return edges_container_.template GetFeatureStorage<Feature>(id);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
template <Component C, typename Feature>
auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                 ExtraStorageT>::GetFeatureExtraStorage() {
  if constexpr (C == Component::Node) {
    return nodes_container_.template GetFeatureExtraStorage<Feature>();
  } else {
    return edges_container_.template GetFeatureExtraStorage<Feature>();
  }
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
template <Component C, typename Feature>
const auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                       ExtraStorageT>::GetFeatureExtraStorage() const {
  if constexpr (C == Component::Node) {
    return nodes_container_.template GetFeatureExtraStorage<Feature>();
  } else {
    return edges_container_.template GetFeatureExtraStorage<Feature>();
  }
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
template <typename Feature>
auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                 ExtraStorageT>::GetFeatureStorage() {
  return features_storage_.template GetFeatureStorage<Feature>();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
template <typename Feature>
const auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                       ExtraStorageT>::GetFeatureStorage() const {
  return features_storage_.template GetFeatureStorage<Feature>();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
template <Component C>
auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                 ExtraStorageT>::GetContainer() {
  if constexpr (C == Component::Node) {
    return nodes_container_;
  } else {
    return edges_container_;
  }
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT>
template <Component C>
const auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                       ExtraStorageT>::GetContainer() const {
  if constexpr (C == Component::Node) {
    return nodes_container_;
  } else {
    return edges_container_;
  }
}

struct DefaultDAGStorage;

template <>
struct LongNameOf<DefaultDAGStorage> {
  using type = DAGStorage<DefaultDAGStorage,
                          ElementsContainer<Component::Node, ElementStorage<Neighbors>>,
                          ElementsContainer<Component::Edge, ElementStorage<Endpoints>>,
                          ExtraStorage<Connections>>;
};

struct DefaultDAGStorage : LongNameOf<DefaultDAGStorage>::type {
  MOVE_ONLY(DefaultDAGStorage);

  using LongNameType = typename LongNameOf<DefaultDAGStorage>::type;
  using LongNameType::LongNameType;

  static inline DefaultDAGStorage EmptyDefault() { return DefaultDAGStorage{}; };
};

#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
           ViewBase>::DAGStorage(NodesContainerT&& nodes_container,
                                 EdgesContainerT&& edges_container,
                                 ExtraStorageT&& features_storage)
    : nodes_container_{std::move(nodes_container)},
      edges_container_{std::move(edges_container)},
      features_storage_{std::move(features_storage)} {}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
template <Component C, typename Feature>
inline constexpr bool DAGStorage<ShortName, NodesContainerT, EdgesContainerT,
                                 ExtraStorageT, ViewBase>::contains_element_feature =
    [] {
      // NOLINTBEGIN
      if constexpr (C == Component::Node) {
        return NodesContainerT::template contains_element_feature<Feature>;
      } else {
        return EdgesContainerT::template contains_element_feature<Feature>;
      }
      // NOLINTEND
    }();

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
typename DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                    ViewBase>::ViewType
DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
           ViewBase>::View() {
  return ViewType{static_cast<Self&>(*this)};
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
typename DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                    ViewBase>::ConstViewType
DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT, ViewBase>::View()
    const {
  return ConstViewType{static_cast<const Self&>(*this)};
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
NodeId DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                  ViewBase>::AppendNode() {
  return nodes_container_.Append();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
EdgeId DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                  ViewBase>::AppendEdge() {
  return edges_container_.Append();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
void DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                ViewBase>::AddNode(NodeId id) {
  nodes_container_.Add(id);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
void DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                ViewBase>::AddEdge(EdgeId id) {
  edges_container_.Add(id);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
size_t DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                  ViewBase>::GetNodesCount() const {
  return nodes_container_.GetCount();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
size_t DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                  ViewBase>::GetEdgesCount() const {
  return edges_container_.GetCount();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
void DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                ViewBase>::InitializeNodes(size_t size) {
  nodes_container_.Initialize(size);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
void DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                ViewBase>::InitializeEdges(size_t size) {
  edges_container_.Initialize(size);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
void DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                ViewBase>::ClearNodes() {
  nodes_container_.Clear();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
void DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                ViewBase>::ClearEdges() {
  edges_container_.Clear();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
template <typename Feature>
auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                 ViewBase>::GetFeatureStorage(NodeId id) {
  return nodes_container_.template GetFeatureStorage<Feature>(id);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
template <typename Feature>
const auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                       ViewBase>::GetFeatureStorage(NodeId id) const {
  return nodes_container_.template GetFeatureStorage<Feature>(id);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
template <typename Feature>
auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                 ViewBase>::GetFeatureStorage(EdgeId id) {
  return edges_container_.template GetFeatureStorage<Feature>(id);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
template <typename Feature>
const auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                       ViewBase>::GetFeatureStorage(EdgeId id) const {
  return edges_container_.template GetFeatureStorage<Feature>(id);
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
template <Component C, typename Feature>
auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                 ViewBase>::GetFeatureExtraStorage() {
  if constexpr (C == Component::Node) {
    return nodes_container_.template GetFeatureExtraStorage<Feature>();
  } else {
    return edges_container_.template GetFeatureExtraStorage<Feature>();
  }
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
template <Component C, typename Feature>
const auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                       ViewBase>::GetFeatureExtraStorage() const {
  if constexpr (C == Component::Node) {
    return nodes_container_.template GetFeatureExtraStorage<Feature>();
  } else {
    return edges_container_.template GetFeatureExtraStorage<Feature>();
  }
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
template <typename Feature>
auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                 ViewBase>::GetFeatureStorage() {
  return features_storage_.template GetFeatureStorage<Feature>();
}

template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
template <typename Feature>
const auto& DAGStorage<ShortName, NodesContainerT, EdgesContainerT, ExtraStorageT,
                       ViewBase>::GetFeatureStorage() const {
  return features_storage_.template GetFeatureStorage<Feature>();
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

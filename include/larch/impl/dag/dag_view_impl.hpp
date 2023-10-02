#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Storage, template <typename, typename> typename Base>
inline constexpr bool DAGView<Storage, Base>::is_mutable = not std::is_const_v<Storage>;

template <typename Storage, template <typename, typename> typename Base>
template <Component C, typename Feature>
inline constexpr bool DAGView<Storage, Base>::contains_element_feature =
    Storage::template contains_element_feature<C, Feature>;

template <typename Storage, template <typename, typename> typename Base>
DAGView<Storage, Base>::DAGView(Storage& dag_storage) : dag_storage_{dag_storage} {}

template <typename Storage, template <typename, typename> typename Base>
DAGView<Storage, Base>::operator DAGView<const Storage, Base>() const {
  return DAGView<const Storage, Base>{dag_storage_};
}

template <typename Storage, template <typename, typename> typename Base>
DAGView<const Storage, Base> DAGView<Storage, Base>::Const() const {
  return DAGView<const Storage, Base>{dag_storage_};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Node, DAGView<Storage, Base>> DAGView<Storage, Base>::Get(
    NodeId id) const {
  return {*this, id};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>> DAGView<Storage, Base>::Get(
    EdgeId id) const {
  return {*this, id};
}

template <typename Storage, template <typename, typename> typename Base>
std::optional<ElementView<Component::Edge, DAGView<Storage, Base>>>
DAGView<Storage, Base>::FindEdge(NodeId parent_id, NodeId child_id) const {
  auto dag = DAGView<const Storage, Base>{dag_storage_};

  auto parent_node = dag.Get(parent_id);
  for (auto edge_id : parent_node.GetChildren()) {
    auto edge = dag.Get(edge_id);
    if (edge.GetChild() == child_id) {
      return {{*this, edge.GetId()}};
    }
  }
  return std::nullopt;
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Node, DAGView<Storage, Base>>
DAGView<Storage, Base>::AppendNode() const {
  NodeId result = dag_storage_.AppendNode();
  return {*this, result};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>>
DAGView<Storage, Base>::AppendEdge() const {
  EdgeId result = dag_storage_.AppendEdge();
  return {*this, result};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Node, DAGView<Storage, Base>> DAGView<Storage, Base>::AddNode(
    NodeId id) {
  dag_storage_.AddNode(id);
  return {*this, id};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>> DAGView<Storage, Base>::AddEdge(
    EdgeId id) {
  dag_storage_.AddEdge(id);
  return {*this, id};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>> DAGView<Storage, Base>::AddEdge(
    EdgeId id, NodeId parent, NodeId child,
    CladeIdx clade) {  // TODO
  dag_storage_.AddEdge(id);
  auto result = Get(id);
  result.Set(parent, child, clade);
  return result;
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>> DAGView<Storage, Base>::AppendEdge(
    NodeId parent, NodeId child,
    CladeIdx clade) const {  // TODO
  auto dag = dag_storage_.View();
  auto result = AppendEdge();
  result.Set(parent, child, clade);
  auto parent_node = dag.Get(parent);
  auto child_node = dag.Get(child);
  parent_node.AddEdge(clade, result.GetId(), true);
  child_node.AddEdge(clade, result.GetId(), false);
  return result;
}

template <typename Storage, template <typename, typename> typename Base>
size_t DAGView<Storage, Base>::GetNodesCount() const {
  return dag_storage_.GetNodesCount();
}

template <typename Storage, template <typename, typename> typename Base>
size_t DAGView<Storage, Base>::GetEdgesCount() const {
  return dag_storage_.GetEdgesCount();
}

template <typename Storage, template <typename, typename> typename Base>
bool DAGView<Storage, Base>::IsEmpty() const {
  return GetNodesCount() == 0 and GetEdgesCount() == 0;
}

template <typename Storage, template <typename, typename> typename Base>
auto DAGView<Storage, Base>::GetNodes() const {
  return ranges::views::indices(GetNodesCount()) |
         ranges::views::transform([*this](size_t i) {
           return ElementView{*this, NodeId{i}};
         });
}

template <typename Storage, template <typename, typename> typename Base>
auto DAGView<Storage, Base>::GetEdges() const {
  return ranges::views::indices(GetEdgesCount()) |
         ranges::views::transform([*this](size_t i) {
           return ElementView{*this, EdgeId{i}};
         });
}

template <typename Storage, template <typename, typename> typename Base>
void DAGView<Storage, Base>::InitializeNodes(size_t size) const {
  dag_storage_.InitializeNodes(size);
}

template <typename Storage, template <typename, typename> typename Base>
void DAGView<Storage, Base>::InitializeEdges(size_t size) const {
  dag_storage_.InitializeEdges(size);
}

template <typename Storage, template <typename, typename> typename Base>
void DAGView<Storage, Base>::ClearNodes() const {
  dag_storage_.ClearNodes();
}

template <typename Storage, template <typename, typename> typename Base>
void DAGView<Storage, Base>::ClearEdges() const {
  dag_storage_.ClearEdges();
}

template <typename Storage, template <typename, typename> typename Base>
template <typename Feature>
auto& DAGView<Storage, Base>::GetFeatureStorage() const {
  return dag_storage_.template GetFeatureStorage<Feature>();
}

template <typename Storage, template <typename, typename> typename Base>
template <typename Feature>
auto& DAGView<Storage, Base>::GetFeatureStorage(NodeId id) const {
  return dag_storage_.template GetFeatureStorage<Feature>(id);
}

template <typename Storage, template <typename, typename> typename Base>
template <typename Feature>
auto& DAGView<Storage, Base>::GetFeatureStorage(EdgeId id) const {
  return dag_storage_.template GetFeatureStorage<Feature>(id);
}

template <typename Storage, template <typename, typename> typename Base>
template <Component C, typename Feature>
auto& DAGView<Storage, Base>::GetFeatureExtraStorage() const {
  return dag_storage_.template GetFeatureExtraStorage<C, Feature>();
}

template <typename Storage, template <typename, typename> typename Base>
const Storage& DAGView<Storage, Base>::GetStorage() const {
  return dag_storage_;
}

template <typename Storage, template <typename, typename> typename Base>
Storage& DAGView<Storage, Base>::GetStorage() {
  return dag_storage_;
}

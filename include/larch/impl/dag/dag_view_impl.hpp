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
DAGView<Storage, Base>::DAGView(Storage& dag_storage)
    : dag_storage_{std::addressof(dag_storage)} {}

template <typename Storage, template <typename, typename> typename Base>
DAGView<Storage, Base>::operator DAGView<const Storage, Base>() const {
  return DAGView<const Storage, Base>{GetStorage()};
}

template <typename Storage, template <typename, typename> typename Base>
DAGView<const Storage, Base> DAGView<Storage, Base>::Const() const {
  return DAGView<const Storage, Base>{GetStorage()};
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
ElementView<Component::Node, DAGView<Storage, Base>>
DAGView<Storage, Base>::AppendNode() const {
  NodeId result = GetStorage().AppendNode();
  return {*this, result};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>>
DAGView<Storage, Base>::AppendEdge() const {
  EdgeId result = GetStorage().AppendEdge();
  return {*this, result};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Node, DAGView<Storage, Base>> DAGView<Storage, Base>::AddNode(
    NodeId id) {
  GetStorage().AddNode(id);
  return {*this, id};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>> DAGView<Storage, Base>::AddEdge(
    EdgeId id) {
  GetStorage().AddEdge(id);
  return {*this, id};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>> DAGView<Storage, Base>::AddEdge(
    EdgeId id, NodeId parent, NodeId child,
    CladeIdx clade) {  // TODO
  GetStorage().AddEdge(id);
  auto result = Get(id);
  result.Set(parent, child, clade);
  return result;
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>> DAGView<Storage, Base>::AppendEdge(
    NodeId parent, NodeId child,
    CladeIdx clade) const {  // TODO
  auto result = AppendEdge();
  result.Set(parent, child, clade);
  return result;
}

template <typename Storage, template <typename, typename> typename Base>
size_t DAGView<Storage, Base>::GetNodesCount() const {
  return GetStorage().GetNodesCount();
}

template <typename Storage, template <typename, typename> typename Base>
size_t DAGView<Storage, Base>::GetEdgesCount() const {
  return GetStorage().GetEdgesCount();
}

template <typename Storage, template <typename, typename> typename Base>
bool DAGView<Storage, Base>::empty() const {
  return GetNodesCount() == 0 and GetEdgesCount() == 0;
}

template <typename Storage, template <typename, typename> typename Base>
auto DAGView<Storage, Base>::GetNodes() const {
  return GetStorage().GetNodes() | Transform::ToNodes(DAGView{*this});
}

template <typename Storage, template <typename, typename> typename Base>
auto DAGView<Storage, Base>::GetEdges() const {
  return GetStorage().GetEdges() | Transform::ToEdges(DAGView{*this});
}

template <typename Storage, template <typename, typename> typename Base>
void DAGView<Storage, Base>::InitializeNodes(size_t size) const {
  GetStorage().InitializeNodes(size);
}

template <typename Storage, template <typename, typename> typename Base>
void DAGView<Storage, Base>::InitializeEdges(size_t size) const {
  GetStorage().InitializeEdges(size);
}

template <typename Storage, template <typename, typename> typename Base>
template <typename Feature>
auto& DAGView<Storage, Base>::GetFeatureStorage() const {
  return GetStorage().template GetFeatureStorage<Feature>();
}

template <typename Storage, template <typename, typename> typename Base>
template <typename Feature>
auto& DAGView<Storage, Base>::GetFeatureStorage(NodeId id) const {
  return GetStorage().template GetFeatureStorage<Feature>(id);
}

template <typename Storage, template <typename, typename> typename Base>
template <typename Feature>
auto& DAGView<Storage, Base>::GetFeatureStorage(EdgeId id) const {
  return GetStorage().template GetFeatureStorage<Feature>(id);
}

template <typename Storage, template <typename, typename> typename Base>
template <Component C, typename Feature>
auto& DAGView<Storage, Base>::GetFeatureExtraStorage() const {
  return GetStorage().template GetFeatureExtraStorage<C, Feature>();
}

template <typename Storage, template <typename, typename> typename Base>
Storage& DAGView<Storage, Base>::GetStorage() const {
  Assert(dag_storage_ != nullptr);
  return *dag_storage_;
}

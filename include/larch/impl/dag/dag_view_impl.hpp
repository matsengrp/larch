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
  // LARCH_DEBUG_USE;
  return DAGView<const Storage, Base>{GetStorage()};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Node, DAGView<Storage, Base>> DAGView<Storage, Base>::Get(
    NodeId id) const {
  // LARCH_DEBUG_USE;
  return {*this, id};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>> DAGView<Storage, Base>::Get(
    EdgeId id) const {
  // LARCH_DEBUG_USE;
  return {*this, id};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Node, DAGView<Storage, Base>>
DAGView<Storage, Base>::AppendNode() const {
  // LARCH_DEBUG_USE;
  NodeId result = GetStorage().AppendNode();
  return {*this, result};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>>
DAGView<Storage, Base>::AppendEdge() const {
  // LARCH_DEBUG_USE;
  EdgeId result = GetStorage().AppendEdge();
  return {*this, result};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Node, DAGView<Storage, Base>> DAGView<Storage, Base>::AddNode(
    NodeId id) {
  // LARCH_DEBUG_USE;
  GetStorage().AddNode(id);
  return {*this, id};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>> DAGView<Storage, Base>::AddEdge(
    EdgeId id) {
  // LARCH_DEBUG_USE;
  GetStorage().AddEdge(id);
  return {*this, id};
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>> DAGView<Storage, Base>::AddEdge(
    EdgeId id, NodeId parent, NodeId child,
    CladeIdx clade) {  // TODO
  // LARCH_DEBUG_USE;
  GetStorage().AddEdge(id);
  auto result = Get(id);
  result.Set(parent, child, clade);
  return result;
}

template <typename Storage, template <typename, typename> typename Base>
ElementView<Component::Edge, DAGView<Storage, Base>> DAGView<Storage, Base>::AppendEdge(
    NodeId parent, NodeId child,
    CladeIdx clade) const {  // TODO
  // LARCH_DEBUG_USE;
  auto result = AppendEdge();
  result.Set(parent, child, clade);
  return result;
}

template <typename Storage, template <typename, typename> typename Base>
size_t DAGView<Storage, Base>::GetNodesCount() const {
  // LARCH_DEBUG_USE;
  return GetStorage().template GetNodesCount<DAGView>();
}

template <typename Storage, template <typename, typename> typename Base>
size_t DAGView<Storage, Base>::GetEdgesCount() const {
  // LARCH_DEBUG_USE;
  return GetStorage().template GetEdgesCount<DAGView>();
}

template <typename Storage, template <typename, typename> typename Base>
bool DAGView<Storage, Base>::empty() const {
  // LARCH_DEBUG_USE;
  return GetNodesCount() == 0 and GetEdgesCount() == 0;
}

template <typename Storage, template <typename, typename> typename Base>
auto DAGView<Storage, Base>::GetNodes() const {
  // LARCH_DEBUG_USE;
  return GetStorage().template GetNodes<DAGView>() | Transform::ToNodes(DAGView{*this});
}

template <typename Storage, template <typename, typename> typename Base>
auto DAGView<Storage, Base>::GetEdges() const {
  // LARCH_DEBUG_USE;
  return GetStorage().template GetEdges<DAGView>() | Transform::ToEdges(DAGView{*this});
}

template <typename Storage, template <typename, typename> typename Base>
void DAGView<Storage, Base>::InitializeNodes(size_t size) const {
  // LARCH_DEBUG_USE;
  GetStorage().InitializeNodes(size);
}

template <typename Storage, template <typename, typename> typename Base>
void DAGView<Storage, Base>::InitializeEdges(size_t size) const {
  // LARCH_DEBUG_USE;
  GetStorage().InitializeEdges(size);
}

template <typename Storage, template <typename, typename> typename Base>
template <typename Feature>
auto DAGView<Storage, Base>::GetFeatureStorage() const {
  // LARCH_DEBUG_USE;
  return GetStorage().template GetFeatureStorage<Feature>();
}

template <typename Storage, template <typename, typename> typename Base>
template <typename Feature>
auto DAGView<Storage, Base>::GetFeatureStorage(NodeId id) const {
  // LARCH_DEBUG_USE;
  return GetStorage().template GetFeatureStorage<Feature>(id);
}

template <typename Storage, template <typename, typename> typename Base>
template <typename Feature>
auto DAGView<Storage, Base>::GetFeatureStorage(EdgeId id) const {
  // LARCH_DEBUG_USE;
  return GetStorage().template GetFeatureStorage<Feature>(id);
}

template <typename Storage, template <typename, typename> typename Base>
template <Component C, typename Feature>
auto DAGView<Storage, Base>::GetFeatureExtraStorage() const {
  // LARCH_DEBUG_USE;
  return GetStorage().template GetFeatureExtraStorage<C, Feature>();
}

template <typename Storage, template <typename, typename> typename Base>
Storage& DAGView<Storage, Base>::GetStorage() const {
  // LARCH_DEBUG_USE;
  Assert(dag_storage_ != nullptr);
  return *dag_storage_;
}

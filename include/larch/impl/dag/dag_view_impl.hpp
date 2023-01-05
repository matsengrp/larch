#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename DagStorageT>
inline constexpr bool DAGView<DagStorageT>::is_mutable =
    not std::is_const_v<DagStorageT>;

template <typename DagStorageT>
template <typename Id, typename Feature>
inline constexpr bool DAGView<DagStorageT>::contains_element_feature =
    DagStorageT::template contains_element_feature<Id, Feature>;

template <typename DagStorageT>
DAGView<DagStorageT>::DAGView(DagStorageT& dag_storage) : dag_storage_{dag_storage} {}

template <typename DagStorageT>
DAGView<DagStorageT>::operator DAGView<const DagStorageT>() const {
  return DAGView<const DagStorageT>{dag_storage_};
}

template <typename DagStorageT>
ElementView<NodeId, DAGView<DagStorageT>> DAGView<DagStorageT>::Get(NodeId id) const {
  return {*this, id};
}

template <typename DagStorageT>
ElementView<EdgeId, DAGView<DagStorageT>> DAGView<DagStorageT>::Get(EdgeId id) const {
  return {*this, id};
}

template <typename DagStorageT>
ElementView<NodeId, DAGView<DagStorageT>> DAGView<DagStorageT>::AppendNode() const {
  NodeId result = dag_storage_.AppendNode();
  return {*this, result};
}

template <typename DagStorageT>
ElementView<EdgeId, DAGView<DagStorageT>> DAGView<DagStorageT>::AppendEdge() const {
  EdgeId result = dag_storage_.AppendEdge();
  return {*this, result};
}

template <typename DagStorageT>
ElementView<NodeId, DAGView<DagStorageT>> DAGView<DagStorageT>::AddNode(NodeId id) {
  dag_storage_.AddNode(id);
  return {*this, id};
}

template <typename DagStorageT>
ElementView<EdgeId, DAGView<DagStorageT>> DAGView<DagStorageT>::AddEdge(
    EdgeId id, NodeId parent, NodeId child,
    CladeIdx clade) {  // TODO
  dag_storage_.AddEdge(id);
  auto result = Get(id);
  result.Set(parent, child, clade);
  return result;
}

template <typename DagStorageT>
ElementView<EdgeId, DAGView<DagStorageT>> DAGView<DagStorageT>::AppendEdge(
    NodeId parent, NodeId child,
    CladeIdx clade) const {  // TODO
  auto result = AppendEdge();
  result.Set(parent, child, clade);
  return result;
}

template <typename DagStorageT>
size_t DAGView<DagStorageT>::GetNodesCount() const {
  return dag_storage_.GetNodesCount();
}

template <typename DagStorageT>
size_t DAGView<DagStorageT>::GetEdgesCount() const {
  return dag_storage_.GetEdgesCount();
}

template <typename DagStorageT>
auto DAGView<DagStorageT>::GetNodes() const {
  return dag_storage_.GetNodes() |
         ranges::views::transform([*this, idx = size_t{}](auto&) mutable {
           return ElementView<NodeId, DAGView<DagStorageT>>{*this, {idx++}};
         });
}

template <typename DagStorageT>
auto DAGView<DagStorageT>::GetEdges() const {
  return dag_storage_.GetEdges() |
         ranges::views::transform([*this, idx = size_t{}](auto&) mutable {
           return ElementView<EdgeId, DAGView<DagStorageT>>{*this, {idx++}};
         });
}

template <typename DagStorageT>
void DAGView<DagStorageT>::InitializeNodes(size_t size) const {
  dag_storage_.InitializeNodes(size);
}

template <typename DagStorageT>
template <typename Feature>
auto& DAGView<DagStorageT>::GetFeatureStorage() const {
  return dag_storage_.template GetFeatureStorage<Feature>();
}

template <typename DagStorageT>
template <typename Feature>
auto& DAGView<DagStorageT>::GetFeatureStorage(NodeId id) const {
  return dag_storage_.template GetFeatureStorage<Feature>(id);
}

template <typename DagStorageT>
template <typename Feature>
auto& DAGView<DagStorageT>::GetFeatureStorage(EdgeId id) const {
  return dag_storage_.template GetFeatureStorage<Feature>(id);
}

template <typename DagStorageT>
template <typename Id, typename Feature>
auto& DAGView<DagStorageT>::GetFeatureExtraStorage() const {
  return dag_storage_.template GetFeatureExtraStorage<Id, Feature>();
}
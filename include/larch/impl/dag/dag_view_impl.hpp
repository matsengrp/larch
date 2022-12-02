#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename DS>
inline constexpr bool DAGView<DS>::is_mutable = not std::is_const_v<DS>;

template <typename DS>
template <typename Id, typename Feature>
inline constexpr bool DAGView<DS>::contains_element_feature =
    DS::template contains_element_feature<Id, Feature>;

template <typename DS>
DAGView<DS>::DAGView(DS& dag_storage) : dag_storage_{dag_storage} {}

template <typename DS>
DAGView<DS>::operator DAGView<const DS>() const {
  return DAGView<const DS>{dag_storage_};
}

template <typename DS>
ElementView<NodeId, DAGView<DS>> DAGView<DS>::Get(NodeId id) const {
  return {*this, id};
}

template <typename DS>
ElementView<EdgeId, DAGView<DS>> DAGView<DS>::Get(EdgeId id) const {
  return {*this, id};
}

template <typename DS>
ElementView<NodeId, DAGView<DS>> DAGView<DS>::AppendNode() const {
  NodeId result = dag_storage_.AppendNode();
  return {*this, result};
}

template <typename DS>
ElementView<EdgeId, DAGView<DS>> DAGView<DS>::AppendEdge() const {
  EdgeId result = dag_storage_.AppendEdge();
  return {*this, result};
}

template <typename DS>
ElementView<NodeId, DAGView<DS>> DAGView<DS>::AddNode(NodeId id) {
  dag_storage_.AddNode(id);
  return {*this, id};
}

template <typename DS>
ElementView<EdgeId, DAGView<DS>> DAGView<DS>::AddEdge(EdgeId id, NodeId parent,
                                                      NodeId child,
                                                      CladeIdx clade) {  // TODO
  dag_storage_.AddEdge(id);
  auto result = Get(id);
  result.Set(parent, child, clade);
  return result;
}

template <typename DS>
ElementView<EdgeId, DAGView<DS>> DAGView<DS>::AppendEdge(
    NodeId parent, NodeId child,
    CladeIdx clade) const {  // TODO
  auto result = AppendEdge();
  result.Set(parent, child, clade);
  return result;
}

template <typename DS>
size_t DAGView<DS>::GetNodesCount() const {
  return dag_storage_.GetNodesCount();
}

template <typename DS>
size_t DAGView<DS>::GetEdgesCount() const {
  return dag_storage_.GetEdgesCount();
}

template <typename DS>
auto DAGView<DS>::GetNodes() const {
  return dag_storage_.GetNodes() |
         ranges::views::transform([*this, idx = size_t{}](auto&) mutable {
           return ElementView<NodeId, DAGView<DS>>{*this, {idx++}};
         });
}

template <typename DS>
auto DAGView<DS>::GetEdges() const {
  return dag_storage_.GetEdges() |
         ranges::views::transform([*this, idx = size_t{}](auto&) mutable {
           return ElementView<EdgeId, DAGView<DS>>{*this, {idx++}};
         });
}

template <typename DS>
void DAGView<DS>::InitializeNodes(size_t size) const {
  dag_storage_.InitializeNodes(size);
}

template <typename DS>
template <typename F>
auto& DAGView<DS>::GetFeatureStorage() const {
  return dag_storage_.template GetFeatureStorage<F>();
}

template <typename DS>
template <typename F>
auto& DAGView<DS>::GetFeatureStorage(NodeId id) const {
  return dag_storage_.template GetFeatureStorage<F>(id);
}

template <typename DS>
template <typename F>
auto& DAGView<DS>::GetFeatureStorage(EdgeId id) const {
  return dag_storage_.template GetFeatureStorage<F>(id);
}

template <typename DS>
template <typename Id, typename F>
auto& DAGView<DS>::GetFeatureExtraStorage() const {
  return dag_storage_.template GetFeatureExtraStorage<Id, F>();
}
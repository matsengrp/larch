#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Storage, typename... Features>
DAGView<Storage, Features...>::DAGView(Storage& storage) : storage_{storage} {}

template <typename Storage, typename... Features>
DAGView<Storage, Features...>::operator Immutable() const {
  return Immutable{storage_};
}

template <typename Storage, typename... Features>
auto DAGView<Storage, Features...>::AddNode(NodeId id) const -> Node {
  Assert(id.value != NoId);
  std::ignore = GetOrInsert(storage_.nodes_.nodes_, id);
  return {*this, id};
}

template <typename Storage, typename... Features>
auto DAGView<Storage, Features...>::AppendNode() const -> Node {
  return AddNode({GetNodesCount()});
}

template <typename Storage, typename... Features>
auto DAGView<Storage, Features...>::AddEdge(EdgeId id, NodeId parent, NodeId child,
                                            CladeIdx clade) const -> Edge {
  Assert(id.value != NoId);
  Assert(parent.value != NoId);
  Assert(child.value != NoId);
  Assert(clade.value != NoId);
  auto& edge_storage = GetOrInsert(storage_.edges_.edges_, id);
  edge_storage.parent_ = parent;
  edge_storage.child_ = child;
  edge_storage.clade_ = clade;
  return {*this, id};
}

template <typename Storage, typename... Features>
auto DAGView<Storage, Features...>::AppendEdge(NodeId parent, NodeId child,
                                               CladeIdx clade) const -> Edge {
  return AddEdge({GetEdgesCount()}, parent, child, clade);
}

// template <typename Storage, typename... Features>
// void DAGView<Storage, Features...>::InitializeNodes(size_t nodes_count) const;

// template <typename Storage, typename... Features>
// void DAGView<Storage, Features...>::BuildConnections() const;

// template <typename Storage, typename... Features>
// void DAGView<Storage, Features...>::BuildConnectionsRaw() const;

template <typename Storage, typename... Features>
auto DAGView<Storage, Features...>::GetNodes() const {
  return storage_.nodes_.nodes_ |
         ranges::views::transform([this, idx = size_t{}](auto&) mutable {
           return Node{*this, {idx++}};
         });
}

template <typename Storage, typename... Features>
auto DAGView<Storage, Features...>::GetEdges() const {
  return storage_.edges_.edges_ |
         ranges::views::transform([this, idx = size_t{}](auto&) mutable {
           return Edge{*this, {idx++}};
         });
}

template <typename Storage, typename... Features>
auto DAGView<Storage, Features...>::Get(NodeId id) const {
  return Node{*this, id};
}

template <typename Storage, typename... Features>
auto DAGView<Storage, Features...>::Get(EdgeId id) const {
  return Edge{*this, id};
}

// template <typename Storage, typename... Features>
// size_t DAGView<Storage, Features...>::GetNodesCount() const;

// template <typename Storage, typename... Features>
// size_t DAGView<Storage, Features...>::GetEdgesCount() const;

// template <typename Storage, typename... Features>
// bool DAGView<Storage, Features...>::IsTree() const;

// template <typename Storage, typename... Features>
// bool DAGView<Storage, Features...>::HaveRoot() const;

template <typename Storage, typename... Features>
auto DAGView<Storage, Features...>::GetRoot() const {
  return Node{*this, {NoId}};
}

// template <typename Storage, typename... Features>
// auto DAGView<Storage, Features...>::GetLeafs() const;

template <typename Storage, typename... Features>
auto& DAGView<Storage, Features...>::GetStorage() const {
  return storage_;
}
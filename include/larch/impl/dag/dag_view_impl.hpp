#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Storage>
DAGView<Storage>::DAGView(Storage& storage) : storage_{storage} {}

template <typename Storage>
DAGView<Storage>::operator Immutable() const {
  return Immutable{storage_};
}

template <typename Storage>
auto DAGView<Storage>::AddNode(NodeId id) -> Node {
  Assert(id.value != NoId);
  std::ignore = GetOrInsert(storage_.nodes_.nodes_, id);
  return {*this, id};
}

template <typename Storage>
auto DAGView<Storage>::AppendNode() -> Node {
  return AddNode({GetNodesCount()});
}

template <typename Storage>
auto DAGView<Storage>::AddEdge(EdgeId id, NodeId parent, NodeId child, CladeIdx clade)
    -> Edge {
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

template <typename Storage>
auto DAGView<Storage>::AppendEdge(NodeId parent, NodeId child, CladeIdx clade) -> Edge {
  return AddEdge({GetEdgesCount()}, parent, child, clade);
}

template <typename Storage>
void DAGView<Storage>::InitializeNodes(size_t nodes_count) {
  storage_.nodes_.nodes_.resize(nodes_count);
}

template <typename Storage>
void DAGView<Storage>::BuildConnections() {
  storage_.root_ = {NoId};
  storage_.leafs_ = {};
  BuildConnectionsRaw();
  for (auto node : GetNodes()) {
    for (auto clade : node.GetClades()) {
      Assert(not clade.empty() && "Empty clade");
    }
    if (node.IsRoot()) {
      Assert(storage_.root_.value == NoId && "Duplicate root");
      storage_.root_ = node;
    }
    if (node.IsLeaf()) {
      storage_.leafs_.push_back(node);
    }
  }
}

template <typename Storage>
void DAGView<Storage>::BuildConnectionsRaw() {
  for (auto& node : storage_.nodes_.nodes_) {
    node.ClearConnections();
  }
  EdgeId edge_id = {0};
  for (auto& edge : storage_.edges_.edges_) {
    Assert(edge.parent_.value != NoId && "Edge has no parent");
    Assert(edge.child_.value != NoId && "Edge has no child");
    Assert(edge.clade_.value != NoId && "Edge has no clade index");
    auto& parent = storage_.nodes_.nodes_.at(edge.parent_.value);
    auto& child = storage_.nodes_.nodes_.at(edge.child_.value);
    parent.AddEdge(edge.clade_, edge_id, true);
    child.AddEdge(edge.clade_, edge_id, false);
    ++edge_id.value;
  }
}

template <typename Storage>
auto DAGView<Storage>::GetNodes() {
  return storage_.nodes_.nodes_ |
         ranges::views::transform([*this, idx = size_t{}](auto&) mutable {
           return Node{*this, {idx++}};
         });
}

template <typename Storage>
auto DAGView<Storage>::GetEdges() {
  return storage_.edges_.edges_ |
         ranges::views::transform([*this, idx = size_t{}](auto&) mutable {
           return Edge{*this, {idx++}};
         });
}

template <typename Storage>
auto DAGView<Storage>::Get(NodeId id) {
  return Node{*this, id};
}

template <typename Storage>
auto DAGView<Storage>::Get(EdgeId id) {
  return Edge{*this, id};
}

template <typename Storage>
size_t DAGView<Storage>::GetNodesCount() {
  return storage_.nodes_.nodes_.size();
}

template <typename Storage>
size_t DAGView<Storage>::GetEdgesCount() {
  return storage_.edges_.edges_.size();
}

template <typename Storage>
bool DAGView<Storage>::IsTree() {
  return GetNodesCount() == GetEdgesCount() + 1;
}

template <typename Storage>
bool DAGView<Storage>::HaveRoot() {
  return storage_.root_.value != NoId;
}

template <typename Storage>
auto DAGView<Storage>::GetRoot() {
  return Node{*this, storage_.root_};
}

template <typename Storage>
auto DAGView<Storage>::GetLeafs() {
  return storage_.nodes_.leafs_ |
         ranges::views::transform([*this, idx = size_t{}](auto&) mutable {
           return Node{*this, {idx++}};
         });
}

template <typename Storage>
auto& DAGView<Storage>::GetStorage() const {
  return storage_;
}

template <typename Storage>
template <typename Feature>
auto& DAGView<Storage>::GetFeatureStorage() const {
  return std::get<Feature>(storage_.features_);
}
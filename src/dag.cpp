#include "dag.hpp"

MutableNode DAG::AddNode(NodeId id) {
  Assert(id.value != NoId);
  [[maybe_unused]] auto& storage = GetOrInsert(nodes_, id);
  return {*this, id};
}

MutableEdge DAG::AddEdge(EdgeId id, NodeId parent, NodeId child, CladeIdx clade) {
  Assert(id.value != NoId);
  Assert(parent.value != NoId);
  Assert(child.value != NoId);
  Assert(clade.value != NoId);
  auto& storage = GetOrInsert(edges_, id);
  storage.Set(parent, child, clade);
  return {*this, id};
}

MutableEdge DAG::AppendEdge(NodeId parent, NodeId child, CladeIdx clade) {
  return AddEdge({GetEdgesCount()}, parent, child, clade);
}

void DAG::InitializeNodes(size_t nodes_count) { nodes_.resize(nodes_count); }

void DAG::BuildConnections() {
  root_ = {NoId};
  leafs_ = {};
  for (auto& node : nodes_) {
    node.ClearConnections();
  }
  EdgeId edge_id = {0};
  for (auto& edge : edges_) {
    Assert(edge.GetParent().value != NoId);
    Assert(edge.GetChild().value != NoId);
    auto& parent = nodes_.at(edge.GetParent().value);
    auto& child = nodes_.at(edge.GetChild().value);
    parent.AddEdge(edge.GetClade(), edge_id, true);
    child.AddEdge(edge.GetClade(), edge_id, false);
    ++edge_id.value;
  }
  for (auto node : GetNodes()) {
    if (node.IsRoot()) {
      root_ = node;
    }
    if (node.IsLeaf()) {
      leafs_.push_back(node);
    }
  }
}

Node DAG::Get(NodeId id) const { return {*this, id}; }
MutableNode DAG::Get(NodeId id) { return {*this, id}; }

Edge DAG::Get(EdgeId id) const { return {*this, id}; }
MutableEdge DAG::Get(EdgeId id) { return {*this, id}; }

size_t DAG::GetNodesCount() const { return nodes_.size(); }
size_t DAG::GetEdgesCount() const { return edges_.size(); }

Node DAG::GetRoot() const { return {*this, root_}; }

MutableNode DAG::GetRoot() { return {*this, root_}; }

std::map<NodeId, NodeId> DAG::ReindexPreOrder() {
  std::map<NodeId, NodeId> index;
  NodeId id = {0};
  auto Reindex = [&](auto& self, Node node) -> void {
    if (index.find(node.GetId()) != index.end()) {
      return;
    }
    index[node.GetId()] = id;
    ++id.value;
    for (Node child : node.GetChildren() | Transform::GetChild()) {
      self(self, child);
    }
  };
  Reindex(Reindex, GetRoot());
  for (auto& edge : edges_) {
    edge.Set(index.at(edge.GetParent()), index.at(edge.GetChild()));
  }
  BuildConnections();
  return index;
}

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
  storage.parent_ = parent;
  storage.child_ = child;
  storage.clade_ = clade;
  return {*this, id};
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
    Assert(edge.parent_.value != NoId);
    Assert(edge.child_.value != NoId);
    auto& parent = nodes_.at(edge.parent_.value);
    auto& child = nodes_.at(edge.child_.value);
    parent.AddEdge(edge.clade_, edge_id, true);
    child.AddEdge(edge.clade_, edge_id, false);
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

Node DAG::GetRoot() const { return {*this, root_}; }

MutableNode DAG::GetRoot() { return {*this, root_}; }

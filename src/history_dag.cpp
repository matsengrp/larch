#include "history_dag.hpp"

MutableNode HistoryDAG::AddNode(NodeId id) {
  assert(id.value != NoId);
  [[maybe_unused]] auto& storage = GetOrInsert(nodes_, id);
  return {*this, id};
}

MutableEdge HistoryDAG::AddEdge(EdgeId id, NodeId parent, NodeId child,
                                CladeIdx clade) {
  assert(id.value != NoId);
  auto& storage = GetOrInsert(edges_, id);
  storage.parent_ = parent;
  storage.child_ = child;
  storage.clade_ = clade;
  return {*this, id};
}

void HistoryDAG::InitializeComponents(size_t nodes_count, size_t edges_count) {
  nodes_.resize(nodes_count);
  edges_.resize(edges_count);
}

void HistoryDAG::BuildConnections() {
  root_ = {NoId};
  leafs_ = {};
  for (auto& node : nodes_) {
    node.ClearConnections();
  }
  EdgeId edge_id = {0};
  for (auto& edge : edges_) {
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

Node HistoryDAG::Get(NodeId id) const { return {*this, id}; }
MutableNode HistoryDAG::Get(NodeId id) { return {*this, id}; }

Edge HistoryDAG::Get(EdgeId id) const { return {*this, id}; }
MutableEdge HistoryDAG::Get(EdgeId id) { return {*this, id}; }

Node HistoryDAG::GetRoot() const { return {*this, root_}; }

MutableNode HistoryDAG::GetRoot() { return {*this, root_}; }

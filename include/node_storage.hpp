#pragma once

class NodeStorage {
  template <typename>
  friend class NodeView;
  friend class DAG;

  void ClearConnections();
  void AddEdge(CladeIdx clade, EdgeId id, bool this_node_is_parent);
  void RemoveEdge(Edge edge, bool this_node_is_parent);

  std::vector<EdgeId> parents_;
  std::vector<std::vector<EdgeId>> clades_;
};

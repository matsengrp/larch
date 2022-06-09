#pragma once

class NodeStorage {
 public:
  const std::vector<EdgeId>& GetParents() const;
  const std::vector<std::vector<EdgeId>>& GetClades() const;

  void ClearConnections();
  void AddEdge(CladeIdx clade, EdgeId id, bool this_node_is_parent);
  void RemoveEdge(Edge edge, bool this_node_is_parent);

 private:
  std::vector<EdgeId> parents_;
  std::vector<std::vector<EdgeId>> clades_;
};

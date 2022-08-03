/**
 * NodeStorage objects keep track of vectors of parent edges, and vectors of
 * vectors of child edges, where two child edges are in the same vector iff
 * they descend from the same child clade.
 *
 * Edges are referenced by their EdgeId, which is their index in DAG.edges_.
*/
#pragma once

class NodeStorage {
 public:
 /**
  * Get vector of parent edges
  */
  const std::vector<EdgeId>& GetParents() const;
  /**
   * Get vectors of child edges corresponding to child clades
   */
  const std::vector<std::vector<EdgeId>>& GetClades() const;

  /**
   * Remove all parent and child edges
   */
  void ClearConnections();
  /**
   * Add a parent or child edge with EdgeId id. If this_node_is_parent is
   * False, clade is ignored.
   */
  void AddEdge(CladeIdx clade, EdgeId id, bool this_node_is_parent);
  /**
   * Remove a parent or child edge.
   */
  void RemoveEdge(Edge edge, bool this_node_is_parent);

 private:
  std::vector<EdgeId> parents_;
  std::vector<std::vector<EdgeId>> clades_;
};

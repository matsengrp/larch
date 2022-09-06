/**
  DAG is the main structure that owns node and edge data, and provides
  various queries.

  Populating with data should be performed by first adding all nodes by the
  AddNode() function, then adding all the edges with AddEdge() and finally
  calling BuildConnections().

  NodeId and EdgeId are strongly typed wrappers around size_t, and data is
  stored internally by the order of its IDs.

  Additional node and edge data may be stored in classes which extend DAG.
  See for example the class MADAG in `mutation_annotated_dag.hpp`.

  GetNodes() and GetEdges() returns a view into the corresponding elements,
  ordered by id.

  TraversePreOrder() returns a view into Nodes in pre-order traversal order.

*/

#pragma once

#include <type_traits>
#include <limits>
#include <iterator>
#include <algorithm>
#include <vector>
#include <string_view>
#include <map>

#include "common.hpp"
#include "node.hpp"
#include "edge.hpp"
#include "node_storage.hpp"
#include "edge_storage.hpp"
#include "traverse_value.hpp"

class DAG {
 public:
  DAG() = default;
  DAG(DAG&&) = default;
  DAG& operator=(DAG&&) = default;
  DAG(const DAG&) = delete;
  DAG& operator=(const DAG&) = delete;

  MutableNode AddNode(NodeId id);

  MutableEdge AddEdge(EdgeId id, NodeId parent, NodeId child, CladeIdx clade);
  MutableEdge AppendEdge(NodeId parent, NodeId child, CladeIdx clade);

  void InitializeNodes(size_t nodes_count);

  /**
   * Properly reference added edges' IDs in node objects, and find root and
   * leaf nodes
   */
  void BuildConnections();

  void BuildConnectionsRaw();

  /**
   * Return a range containing NodeId's for each node in the DAG
   * @{
   */
  inline auto GetNodes() const;
  inline auto GetNodes();
  /** @} */

  /**
   * Return a range containing EdgeId's for each edge in the DAG
   * @{
   */
  inline auto GetEdges() const;
  inline auto GetEdges();
  /** @} */

  /**
   * Get a Node object by its NodeId
   * @{
   */
  Node Get(NodeId id) const;
  MutableNode Get(NodeId id);
  /** @} */
  /**
   * Get an Edge object by its EdgeId
   * @{
   */
  Edge Get(EdgeId id) const;
  MutableEdge Get(EdgeId id);
  /** @} */

  void DAG::DepthFirstExpansionHelper(NodeId node, std::vector<NodeId>& vec) const;
  std::vector<NodeId> DAG::DepthFirstExpansion(NodeId node) const {

  size_t GetNodesCount() const;
  size_t GetEdgesCount() const;

  bool IsTree() const;

  bool HaveRoot() const;
  Node GetRoot() const;
  MutableNode GetRoot();

  /**
   * Return a range containing leaf Nodes in the DAG
   * @{
   */
  inline auto GetLeafs() const;
  inline auto GetLeafs();
  /** @} */

  /**
   * Return a range containing a preordering of Nodes in the DAG
   * @{
   */
  inline auto TraversePreOrder() const;
  inline auto TraversePreOrder();
  /** @} */

  /**
   * Return a range containing a postordering of Nodes in the DAG
   * @{
   */
  inline auto TraversePostOrder() const;
  inline auto TraversePostOrder();
  /** @} */

  /**
   * Change node IDs so that they are pre-ordered, and return a
   * map from old NodeIds to new NodeIds.
   */
  std::map<NodeId, NodeId> ReindexPreOrder();

 private:
  template <typename>
  friend class NodeView;
  template <typename>
  friend class EdgeView;

  std::vector<NodeStorage> nodes_;
  std::vector<EdgeStorage> edges_;

  NodeId root_ = {NoId};
  std::vector<NodeId> leafs_;
};

#include "pre_order_iterator.hpp"
#include "post_order_iterator.hpp"

#include "impl/node_impl.hpp"
#include "impl/dag_impl.hpp"
#include "impl/traverse_value_impl.hpp"

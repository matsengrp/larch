/*
  HistoryDAG is the main structure that owns node and edge data, and provides
  various queries.

  Populating with data should be performed by first adding all nodes by the
  AddNode() function, then adding all the edges with AddEdge() and finally
  calling BuildConnections().

  NodeId and EdgeId are strongly typed wrappers around size_t, and data is
  stored internally by the order of its IDs.

  GetNodes() and GetEdges() returns a view into the corresponding elements,
  ordered by id.

  TraversePreOrder() returns a view into Nodes in pre-order traversal order.

*/
#pragma once

#include <type_traits>
#include <limits>
#include <iterator>
#include <algorithm>
#include <cassert>
#include <vector>
#include <string_view>

#include "history_dag_node_storage.hpp"
#include "history_dag_edge_storage.hpp"
#include "traverse_value.hpp"
#include "counter_map.hpp"

class HistoryDAG {
 public:
  HistoryDAG() = default;
  HistoryDAG(HistoryDAG&&) = default;
  HistoryDAG& operator=(HistoryDAG&&) = default;

  using Weight = double;
  using ArbitraryPrecisionInteger = long;

  MutableNode AddNode(NodeId id);

  MutableEdge AddEdge(EdgeId id, Node parent, Node child, CladeIdx clade);
  MutableEdge AddEdge(EdgeId id, NodeId parent, NodeId child, CladeIdx clade);

  void BuildConnections();

  inline auto GetNodes() const;
  inline auto GetNodes();
  inline auto GetEdges() const;
  inline auto GetEdges();

  Node GetNode(NodeId id) const;
  MutableNode GetNode(NodeId id);
  Edge GetEdge(EdgeId id) const;
  MutableEdge GetEdge(EdgeId id);

  Node GetRoot() const;
  MutableNode GetRoot();

  inline auto GetLeafs() const;
  inline auto GetLeafs();

  inline auto TraversePreOrder() const;
  inline auto TraversePreOrder();
  inline auto TraversePostOrder() const;
  inline auto TraversePostOrder();

  ArbitraryPrecisionInteger CountHistories() const;
  void WriteProtobuf(std::string_view filename) const;
  HistoryDAG SampleHistory() const;
  HistoryDAG FindHistoryByIndex(ArbitraryPrecisionInteger) const;
  bool IsCladeTree() const;
  void AddAllAllowedEdges();

 private:
  template <typename>
  friend class NodeView;
  template <typename>
  friend class EdgeView;

  std::vector<NodeStorage> nodes_;

  std::vector<EdgeStorage<Weight>> edges_;

  NodeId root_ = {NoId};
  std::vector<NodeId> leafs_;
};

#include "history_dag_node.hpp"
#include "history_dag_edge.hpp"
#include "pre_order_iterator.hpp"
#include "post_order_iterator.hpp"
#include "impl/history_dag_node_impl.hpp"
#include "impl/history_dag_edge_impl.hpp"
#include "impl/history_dag_impl.hpp"
#include "impl/pre_order_iterator_impl.hpp"
#include "impl/post_order_iterator_impl.hpp"
#include "impl/traverse_value_impl.hpp"

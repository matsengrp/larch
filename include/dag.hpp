/*
  DAG is the main structure that owns node and edge data, and provides
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
#include <vector>
#include <string_view>

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

  using Weight = double;
  using ArbitraryPrecisionInteger = long;

  MutableNode AddNode(NodeId id);

  MutableEdge AddEdge(EdgeId id, NodeId parent, NodeId child, CladeIdx clade);

  void InitializeNodes(size_t nodes_count);

  void BuildConnections();

  inline auto GetNodes() const;
  inline auto GetNodes();
  inline auto GetEdges() const;
  inline auto GetEdges();

  Node Get(NodeId id) const;
  MutableNode Get(NodeId id);
  Edge Get(EdgeId id) const;
  MutableEdge Get(EdgeId id);

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
  DAG SampleHistory() const;
  DAG FindHistoryByIndex(ArbitraryPrecisionInteger) const;
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

#include "pre_order_iterator.hpp"
#include "post_order_iterator.hpp"
#include "impl/node_impl.hpp"
#include "impl/edge_impl.hpp"
#include "impl/dag_impl.hpp"
#include "impl/pre_order_iterator_impl.hpp"
#include "impl/post_order_iterator_impl.hpp"
#include "impl/traverse_value_impl.hpp"

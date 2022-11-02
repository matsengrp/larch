#pragma once

#include "larch/common.hpp"
#include <optional>
struct NodeId {
  size_t value = NoId;
};

inline bool operator==(NodeId lhs, NodeId rhs) { return lhs.value == rhs.value; }

inline bool operator<(NodeId lhs, NodeId rhs) { return lhs.value < rhs.value; }

template <>
struct std::hash<NodeId> {
  size_t operator()(NodeId id) const noexcept { return id.value; }
};

/*
  NodeView is a lightweight view object into the internal node storage of a
  DAG. It is meant to be cheaply passed by value, and behaves as a
  reference into raw storage that conveniently enriches its API.
*/
template <typename T>
class NodeView {
 public:
  constexpr static const bool is_mutable = std::is_same_v<T, DAG&>;
  using NodeType = std::conditional_t<is_mutable, MutableNode, Node>;
  using EdgeType = std::conditional_t<is_mutable, MutableEdge, Edge>;
  NodeView(T dag, NodeId id);
  operator Node() const;
  operator NodeId() const;
  /**
   * Return DAG-like object containing this node
   */
  T GetDAG() const;
  NodeId GetId() const;
  /**
   * Return a range containing parent Edge objects
   */
  auto GetParents() const;
  /**
   * Return a range containing clades. A clade is a range containing child
   * Edges
   */
  auto GetClades() const;
  /**
   * Return a clade by it's index. A clade is a range containing child Edges
   */
  auto GetClade(CladeIdx clade) const;
  /**
   * Return the count of child clades
   */
  size_t GetCladesCount() const;
  /**
   * Return a range containing child Edges
   */
  auto GetChildren() const;
  /**
   * Return a single parent edge of this node
   */
  EdgeType GetSingleParent() const;

  EdgeType GetFirstChild() const;
  auto GetFirstClade() const;

  /**
   * Checks if node has no parents
   */
  bool IsRoot() const;
  /**
   * Checks if node has no children
   */
  bool IsLeaf() const;
  void AddParentEdge(Edge edge) const;
  void AddChildEdge(Edge edge) const;
  void RemoveParentEdge(Edge edge) const;

  const std::optional<std::string>& GetSampleId() const;
  void SetSampleId(std::optional<std::string>&& sample_id) const;

 private:
  auto& GetStorage() const;
  T dag_;
  const NodeId id_;
};

template <typename T>
inline bool operator==(NodeView<T> lhs, NodeView<T> rhs) {
  return std::addressof(lhs.GetDAG()) == std::addressof(rhs.GetDAG) &&
         lhs.GetId() == rhs.GetId();
}

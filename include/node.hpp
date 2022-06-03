/*
  Node is a lightweight view object into the internal node storage of a
  DAG. It is meant to be cheaply passed by value, and behaves as a
  reference into raw storage that conveniently enriches it's API.
*/
#pragma once

struct NodeId {
  size_t value = NoId;
};

inline bool operator==(NodeId lhs, NodeId rhs) { return lhs.value == rhs.value; }

inline bool operator<(NodeId lhs, NodeId rhs) { return lhs.value < rhs.value; }

template <>
struct std::hash<NodeId> {
  size_t operator()(NodeId id) const noexcept { return id.value; }
};

template <typename T>
class NodeView {
 public:
  constexpr static const bool is_mutable = std::is_same_v<T, DAG&>;
  using NodeType = std::conditional_t<is_mutable, MutableNode, Node>;
  using EdgeType = std::conditional_t<is_mutable, MutableEdge, Edge>;
  NodeView(T dag, NodeId id);
  operator Node() const;
  operator NodeId() const;
  T GetDAG() const;
  NodeId GetId() const;
  auto GetParents() const;
  auto GetClades() const;
  auto GetChildren() const;
  EdgeType GetSingleParent() const;
  bool IsRoot() const;
  bool IsLeaf() const;
  void AddParentEdge(Edge edge);
  void AddChildEdge(Edge edge);
  void RemoveParentEdge(Edge edge);

 private:
  template <typename U>
  friend bool operator==(NodeView<U>, NodeView<U>);
  template <typename>
  friend class NodeView;
  auto& GetStorage() const;
  T dag_;
  const NodeId id_;
};

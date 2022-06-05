/*
  Edge is a lightweight view object into the internal edge storage of a
  DAG. It is meant to be cheaply passed by value, and behaves as a
  reference into raw storage that conveniently enriches it's API.
*/
#pragma once

#include <optional>

#include "common.hpp"

struct EdgeId {
  size_t value = NoId;
};

inline bool operator==(EdgeId lhs, EdgeId rhs) { return lhs.value == rhs.value; }

inline bool operator<(EdgeId lhs, EdgeId rhs) { return lhs.value < rhs.value; }

struct CladeIdx {
  size_t value = NoId;
};

inline bool operator==(CladeIdx lhs, CladeIdx rhs) { return lhs.value == rhs.value; }

inline bool operator<(CladeIdx lhs, CladeIdx rhs) { return lhs.value < rhs.value; }

template <typename T>
class EdgeView {
 public:
  constexpr static const bool is_mutable = std::is_same_v<T, DAG&>;
  using NodeType = std::conditional_t<is_mutable, MutableNode, Node>;
  EdgeView(T dag, EdgeId id);
  operator Edge() const;
  operator EdgeId() const;
  operator CladeIdx() const;
  T GetDAG() const;
  EdgeId GetId() const;
  NodeType GetParent() const;
  NodeType GetChild() const;
  CladeIdx GetClade() const;
  NodeId GetParentId() const;
  NodeId GetChildId() const;
  std::pair<NodeId, NodeId> GetNodeIds() const;
  bool IsRoot() const;
  bool IsLeaf() const;
  double GetProbability() const;
  const auto& GetWeight() const;
  std::optional<EdgeView> FindNextSibling() const;

 private:
  template <typename U>
  friend bool operator==(EdgeView<U>, EdgeView<U>);
  const auto& GetStorage() const;
  auto& GetStorage();

  T dag_;
  const EdgeId id_;
};

template <typename T>
inline bool operator==(EdgeView<T> lhs, EdgeView<T> rhs) {
  return std::addressof(lhs.dag_) == std::addressof(rhs.dag_) && lhs.id_ == rhs.id_;
}

namespace std {
template <typename T>
struct tuple_size<::EdgeView<T>> : integral_constant<size_t, 2> {};

template <size_t Index, typename T>
struct tuple_element<Index, ::EdgeView<T>>
    : tuple_element<Index, tuple<NodeView<T>, NodeView<T>>> {};
}  // namespace std

template <std::size_t Index, typename T>
std::tuple_element_t<Index, EdgeView<T>> get(EdgeView<T> edge) {
  if constexpr (Index == 0) return edge.GetParent();
  if constexpr (Index == 1) return edge.GetChild();
}

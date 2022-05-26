/*
  Edge is a lightweight view object into the internal edge storage of a
  HistoryDAG. It is meant to be cheaply passed by value, and behaves as a
  reference into raw storage that conveniently enriches it's API.
*/
#pragma once

#include <optional>

#include "history_dag_common.hpp"

class HistoryDAG;

template <typename T>
class EdgeView {
 public:
  constexpr static const bool is_mutable = std::is_same_v<T, HistoryDAG&>;
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

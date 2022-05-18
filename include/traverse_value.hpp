#pragma once

#include "history_dag_common.hpp"

template <typename T>
class TraverseValue {
 public:
  TraverseValue(T dag, NodeId node, EdgeId edge);

  NodeView<T> GetNode() const;
  EdgeView<T> GetEdge() const;

  operator MutableNode() const;
  operator MutableEdge() const;
  operator Node() const;
  operator Edge() const;

  bool IsRoot() const;
  bool IsLeaf() const;

 private:
  T dag_;
  const NodeId node_;
  const EdgeId edge_;
};

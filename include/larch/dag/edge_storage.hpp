#pragma once

#include "larch/dag/node.hpp"
#include "larch/dag/edge.hpp"

class EdgeStorage {
 public:
  NodeId GetParent() const;
  NodeId GetChild() const;
  CladeIdx GetClade() const;
  void Set(NodeId parent, NodeId child, CladeIdx clade);
  void Set(NodeId parent, NodeId child);

 private:
  NodeId parent_;
  NodeId child_;
  CladeIdx clade_;
};

#pragma once

#include "larch/dag/dag.hpp"

struct LCA {
  NodeId lca;
  std::vector<EdgeId> path0;
  std::vector<EdgeId> path1;
};

template <typename Node>
LCA FindLCA(Node n0, Node n1);

#include "larch/impl/spr/lca_impl.hpp"
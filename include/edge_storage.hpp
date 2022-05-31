#pragma once

template <typename Weight>
class EdgeStorage {
  template <typename>
  friend class EdgeView;
  friend class DAG;

  NodeId parent_;
  NodeId child_;
  CladeIdx clade_;
  double probability_ = 0.0;
  Weight weight_ = {};
};

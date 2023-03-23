#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Basic per-edge feature.
 */
struct Endpoints {
  NodeId parent_;
  NodeId child_;
  CladeIdx clade_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<Endpoints, CRTP, Tag> {
  auto GetParent() const;
  auto GetChild() const;
  CladeIdx GetClade() const;
  NodeId GetParentId() const;
  NodeId GetChildId() const;
  std::pair<NodeId, NodeId> GetNodeIds() const;
  bool IsRoot() const;
  bool IsLeaf() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<Endpoints, CRTP, Tag> {
  void Set(NodeId parent, NodeId child, CladeIdx clade) const;
  void SetClade(CladeIdx clade) const;
};

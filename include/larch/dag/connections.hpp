#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

#include "larch/contiguous_set.hpp"

/**
 * Basic per-DAG feature.
 */
struct Connections {
  MOVE_ONLY_DEF_CTOR(Connections);
  explicit Connections(NodeId root_node_id) : root_{root_node_id} {}
  NodeId root_;
  std::vector<NodeId> leafs_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<Connections, CRTP, Tag> {
  bool IsTree() const;
  bool HaveRoot() const;
  auto GetRoot() const;

  /**
   * Return a range containing leaf Nodes in the DAG
   */
  auto GetLeafs() const;
  auto GetLeafsCount() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<Connections, CRTP, Tag> {
  /**
   * Properly reference added edges' IDs in node objects, and find root and
   * leaf nodes
   */
  void BuildConnections() const;
  void BuildConnectionsRaw() const;
  void AddLeaf(NodeId id) const;

  ContiguousMap<ContiguousSet<NodeId>, ContiguousSet<NodeId>> BuildCladeUnionMap()
      const;
  void MakeComplete() const;
  void ClearConnections() const;
};

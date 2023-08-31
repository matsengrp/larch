#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Basic per-DAG feature.
 */
struct Connections {
  NodeId root_ = {NoId};
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
};

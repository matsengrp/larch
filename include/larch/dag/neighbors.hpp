#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Basic per-node feature.
 */
struct Neighbors {
  std::vector<EdgeId> parents_;
  std::vector<std::vector<EdgeId>> clades_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<Neighbors, CRTP, Tag> {
  auto GetParents() const;
  auto GetClades() const;
  auto GetClade(CladeIdx clade) const;
  size_t GetParentsCount() const;
  size_t GetCladesCount() const;
  auto GetChildren() const;
  auto GetSingleParent() const;
  auto GetFirstParent() const;
  auto GetFirstChild() const;
  auto GetFirstClade() const;
  bool IsRoot() const;
  bool IsLeaf() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<Neighbors, CRTP, Tag> {
  void ClearConnections() const;
  void AddEdge(CladeIdx clade, EdgeId id, bool this_node_is_parent) const;
};
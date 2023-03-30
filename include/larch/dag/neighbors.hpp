#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Basic per-node feature.
 */
struct Neighbors {
  std::vector<EdgeId> parents_;
  std::vector<std::vector<EdgeId>> clades_;
  std::vector<std::vector<NodeId>> leafs_below_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<Neighbors, CRTP, Tag> {
  auto GetParents() const;
  auto GetClades() const;
  auto GetClade(CladeIdx clade) const;
  size_t GetParentsCount() const;
  size_t GetCladesCount() const;
  auto GetChildren() const;
  size_t GetChildrenCount() const;
  auto GetSingleParent() const;
  auto GetFirstParent() const;
  auto GetFirstChild() const;
  auto GetFirstClade() const;
  bool IsRoot() const;
  bool IsLeaf() const;
  auto GetParentNodes() const;
  auto GetChildNodes() const;
  bool ContainsParent(NodeId node) const;
  bool ContainsChild(NodeId node) const;
  std::string ParentsToString() const;
  std::string ChildrenToString() const;
  auto GetLeafsBelow() const;
  void Validate(bool recursive = false) const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<Neighbors, CRTP, Tag> {
  void ClearConnections() const;
  void AddEdge(CladeIdx clade, EdgeId id, bool this_node_is_parent) const;
  void RemoveParent(EdgeId edge) const;
  void ChangeParent(EdgeId from, EdgeId to) const;
  void SetSingleParent(EdgeId parent) const;
  void RemoveChild(CladeIdx clade, EdgeId child) const;
  void ChangeChild(CladeIdx clade, EdgeId from, EdgeId to) const;
  void CalculateLeafsBelow() const;
};

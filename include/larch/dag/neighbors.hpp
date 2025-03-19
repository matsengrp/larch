#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Basic per-node feature.
 */
struct Neighbors {
  MOVE_ONLY(Neighbors);
  Neighbors() = default;
};

struct DAGNeighbors : Neighbors {
  MOVE_ONLY(DAGNeighbors);
  DAGNeighbors() = default;

  inline DAGNeighbors Copy() const {
    DAGNeighbors result;
    result.parents_ = parents_;
    result.clades_ = clades_;
    result.leafs_below_ = leafs_below_;
    return result;
  }

  template <typename CRTP>
  auto GetParents(const CRTP*) const {
    return parents_ | ranges::view::all;
  }
  template <typename CRTP>
  auto GetClades(const CRTP*) const {
    return clades_ | ranges::view::all;
  }
  template <typename CRTP>
  auto GetLeafsBelow(const CRTP*) const {
    return leafs_below_ | ranges::view::all;
  }

  template <typename CRTP>
  auto& GetParentsMutable(const CRTP*) {
    return parents_;
  }
  template <typename CRTP>
  auto& GetCladesMutable(const CRTP*) {
    return clades_;
  }
  template <typename CRTP>
  auto& GetLeafsBelowMutable(const CRTP*) {
    return leafs_below_;
  }

 private:
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
  auto GetSingleParent() const;
  auto GetFirstParent() const;
  auto GetFirstChild() const;
  auto GetFirstClade() const;
  bool IsUA() const;
  bool IsTreeRoot() const;
  bool IsLeaf() const;
  auto GetLeafsBelow() const;
  void Validate(bool recursive = false, bool allow_dag = false) const;
  auto GetParentNodes() const;
  auto GetChildNodes() const;
  bool ContainsParent(NodeId node) const;
  bool ContainsChild(NodeId node) const;
  std::string ParentsToString() const;
  std::string ChildrenToString() const;
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

template <typename CRTP, typename Tag>
struct FeatureConstView<DAGNeighbors, CRTP, Tag>
    : FeatureConstView<Neighbors, CRTP, Tag> {};

template <typename CRTP, typename Tag>
struct FeatureMutableView<DAGNeighbors, CRTP, Tag>
    : FeatureMutableView<Neighbors, CRTP, Tag> {};

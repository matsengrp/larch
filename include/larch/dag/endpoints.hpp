#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Basic per-edge feature.
 */
struct Endpoints {
  MOVE_ONLY(Endpoints);
  Endpoints() = default;
};

struct DAGEndpoints : Endpoints {
  inline DAGEndpoints Copy() const {
    DAGEndpoints result;
    result.parent_ = parent_;
    result.child_ = child_;
    result.clade_ = clade_;
    return result;
  }

  template <typename CRTP>
  NodeId GetParent(const CRTP*) const {
    return parent_;
  }
  template <typename CRTP>
  NodeId GetChild(const CRTP*) const {
    return child_;
  }
  template <typename CRTP>
  CladeIdx GetClade(const CRTP*) const {
    return clade_;
  }

  template <typename CRTP>
  void SetParent(const CRTP*, NodeId parent) {
    parent_ = parent;
  }
  template <typename CRTP>
  void SetChild(const CRTP*, NodeId child) {
    child_ = child;
  }
  template <typename CRTP>
  void SetClade(const CRTP*, CladeIdx clade) {
    clade_ = clade;
  }

 private:
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
  bool IsUA() const;
  bool IsTreeRoot() const;
  bool IsLeaf() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<Endpoints, CRTP, Tag> {
  void Set(NodeId parent, NodeId child, CladeIdx clade) const;
  void SetParent(NodeId parent) const;
  void SetChild(NodeId child) const;
  void SetClade(CladeIdx clade) const;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<DAGEndpoints, CRTP, Tag>
    : FeatureConstView<Endpoints, CRTP, Tag> {};

template <typename CRTP, typename Tag>
struct FeatureMutableView<DAGEndpoints, CRTP, Tag>
    : FeatureMutableView<Endpoints, CRTP, Tag> {};

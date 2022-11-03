#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename DAGType, typename... Features>
class NodeView
    : public std::conditional_t<
          DAGType::is_mutable, FeatureWriter<Features, NodeView<DAGType, Features...>>,
          FeatureReader<Features, NodeView<DAGType, Features...>>>... {
 public:
  constexpr static const bool is_mutable = DAGType::is_mutable;
  NodeView(const DAGType& dag, NodeId id);
  operator NodeView<typename DAGType::Immutable, Features...>() const;
  operator NodeId() const;
  auto& GetDAG() const;
  NodeId GetId() const;
  auto GetParents() const;
  auto GetClades() const;
  auto GetClade(CladeIdx clade) const;
  size_t GetCladesCount() const;
  auto GetChildren() const;
  auto GetSingleParent() const;
  auto GetFirstChild() const;
  auto GetFirstClade() const;
  bool IsRoot() const;
  bool IsLeaf() const;
  template <typename EdgeView>
  void AddParentEdge(EdgeView edge) const;
  template <typename EdgeView>
  void AddChildEdge(EdgeView edge) const;
  template <typename EdgeView>
  void RemoveParentEdge(EdgeView edge) const;

 private:
  DAG_FEATURE_FRIENDS;
  auto& GetStorage() const;
  const DAGType& dag_;
  NodeId id_;
};
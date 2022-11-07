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
  NodeView(DAGType dag, NodeId id);
  operator NodeView<typename DAGType::Immutable, Features...>();
  operator NodeId();
  auto& GetDAG();
  NodeId GetId();
  auto GetParents();
  auto GetClades();
  auto GetClade(CladeIdx clade);
  size_t GetCladesCount();
  auto GetChildren();
  auto GetSingleParent();
  auto GetFirstChild();
  auto GetFirstClade();
  bool IsRoot();
  bool IsLeaf();
  template <typename EdgeView>
  void AddParentEdge(EdgeView edge);
  template <typename EdgeView>
  void AddChildEdge(EdgeView edge);
  template <typename EdgeView>
  void RemoveParentEdge(EdgeView edge);

 private:
  DAG_FEATURE_FRIENDS;
  auto& GetStorage() const;
  DAGType dag_;
  NodeId id_;
};
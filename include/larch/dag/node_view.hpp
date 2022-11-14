#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename DAG>
class NodeView
    : public DAG::StorageType::NodesContainerType::StorageType::template ViewBase<DAG>,
      public DAG::StorageType::NodesContainerType::template ViewBase<DAG> {
 public:
  constexpr static const bool is_mutable = DAG::is_mutable;
  NodeView(DAG dag, NodeId id);
  operator NodeView<typename DAG::Immutable>();
  operator NodeId();
  template <typename Feature>
  auto& Get();
  template <typename Feature>
  void Set(Feature&& feature);
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
  DAG dag_;
  NodeId id_;
};
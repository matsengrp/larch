#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename DAG>
class Fragment {
 public:
  MOVE_ONLY(Fragment);
  template <Component C, typename Feature>
  static const bool contains_element_feature =
      DAG::template contains_element_feature<C, Feature>;

  Fragment(DAG dag, std::vector<NodeId>&& nodes, std::vector<EdgeId>&& edges, NodeId root_node_id);

  void AssertUA() const;
  size_t GetNodesCount() const;
  size_t GetEdgesCount() const;
  auto Get(NodeId id) const;
  auto Get(EdgeId id) const;
  auto GetNodes() const;
  auto GetEdges() const;
  auto GetRoot() const;
  bool IsTree() const;

 private:
  DAG dag_;
  const std::vector<NodeId> nodes_;
  const std::vector<EdgeId> edges_;
  const NodeId root_node_id_;
};

#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename DAG>
class Fragment {
 public:
  MOVE_ONLY(Fragment);
  template <typename Id, typename Feature>
  static const bool contains_element_feature =
      DAG::template contains_element_feature<Id, Feature>;

  Fragment(DAG dag, std::vector<NodeId>&& nodes, std::vector<EdgeId>&& edges);

  void AssertUA() const;
  size_t GetNodesCount() const;
  size_t GetEdgesCount() const;
  auto Get(NodeId id) const;
  auto Get(EdgeId id) const;
  auto GetNodes() const;
  auto GetEdges() const;
  auto GetRoot() const;

 private:
  DAG dag_;
  const std::vector<NodeId> nodes_;
  const std::vector<EdgeId> edges_;
};

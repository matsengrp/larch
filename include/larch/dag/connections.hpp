#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

struct Connections {
  NodeId root_ = {NoId};
  std::vector<NodeId> leafs_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<Connections, CRTP, Tag> {
  bool IsTree() const;
  bool HaveRoot() const;
  auto GetRoot() const;
  auto GetLeafs() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<Connections, CRTP, Tag> {
  void BuildConnections() const;
  void BuildConnectionsRaw() const;
};
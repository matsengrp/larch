#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

struct MappedNodes {
  NodeId original_id_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<MappedNodes, CRTP, Tag> {
  NodeId GetOriginalId() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<MappedNodes, CRTP, Tag> {
  void SetOriginalId(NodeId id) const;
};

template <>
struct ExtraFeatureStorage<MappedNodes> {
  ContiguousMap<NodeId, NodeId> reverse_map_;
};

template <typename CRTP>
struct ExtraFeatureConstView<MappedNodes, CRTP> {
  auto GetMappedNode(NodeId original_id) const;
};

template <typename CRTP>
struct ExtraFeatureMutableView<MappedNodes, CRTP> {
  auto GetMutableMappedNode(NodeId original_id) const;
};

template <typename DAG>
struct MappedNodesStorage
    : ExtendDAGStorage<MappedNodesStorage<DAG>, DAG, Extend::Nodes<MappedNodes>> {
  using ExtendDAGStorage<MappedNodesStorage<DAG>, DAG,
                         Extend::Nodes<MappedNodes>>::ExtendDAGStorage;
};

template <typename DAG>
MappedNodesStorage<DAG> AddMappedNodes(DAG&& dag) {
  return MappedNodesStorage<DAG>{std::forward<DAG>(dag)};
}

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
    : ExtendStorageType<MappedNodesStorage<DAG>, DAG, Extend::Nodes<MappedNodes>> {
  using ExtendStorageType<MappedNodesStorage<DAG>, DAG,
                          Extend::Nodes<MappedNodes>>::ExtendDAGStorage;
  MappedNodesStorage(MappedNodesStorage&& other) = default;
  MappedNodesStorage& operator=(MappedNodesStorage&& other) {
    static_cast<
        ExtendStorageType<MappedNodesStorage<DAG>, DAG, Extend::Nodes<MappedNodes>>&>(
        *this) = std::move(other);
    return *this;
  }
};

template <typename DAG, typename = std::enable_if_t<DAG::role == Role::Storage>>
MappedNodesStorage<DAG> AddMappedNodes(DAG&& dag) {
  return MappedNodesStorage<DAG>::Consume(std::move(dag));
}

template <typename DAG, typename = std::enable_if_t<DAG::role == Role::View>>
MappedNodesStorage<DAG> AddMappedNodes(const DAG& dag) {
  return MappedNodesStorage<DAG>::FromView(std::move(dag));
}

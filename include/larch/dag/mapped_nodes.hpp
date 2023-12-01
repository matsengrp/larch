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
  MOVE_ONLY_DEF_CTOR(ExtraFeatureStorage);
  IdContainer<NodeId, NodeId, IdContinuity::Sparse, Ordering::Ordered> reverse_map_;
};

template <typename CRTP>
struct ExtraFeatureConstView<MappedNodes, CRTP> {
  auto GetMappedNode(NodeId original_id) const;
};

template <typename CRTP>
struct ExtraFeatureMutableView<MappedNodes, CRTP> {
  auto GetMutableMappedNode(NodeId original_id) const;
};

template <typename Target>
struct MappedNodesStorage;

template <typename Target>
struct LongNameOf<MappedNodesStorage<Target>> {
  using type =
      ExtendStorageType<MappedNodesStorage<Target>, Target, Extend::Nodes<MappedNodes>>;
};

template <typename Target>
struct MappedNodesStorage : LongNameOf<MappedNodesStorage<Target>>::type {
  SHORT_NAME(MappedNodesStorage);
};

template <typename DAG, typename = std::enable_if_t<DAG::role == Role::Storage>>
MappedNodesStorage<DAG> AddMappedNodes(DAG&& dag) {
  return MappedNodesStorage<DAG>::Consume(std::move(dag));
}

template <typename DAG, typename = std::enable_if_t<DAG::role == Role::View>>
MappedNodesStorage<DAG> AddMappedNodes(const DAG& dag) {
  return MappedNodesStorage<DAG>::FromView(dag);
}

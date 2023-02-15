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
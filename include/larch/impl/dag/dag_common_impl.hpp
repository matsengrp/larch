#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

bool operator==(NodeId lhs, NodeId rhs) { return lhs.value == rhs.value; }

bool operator<(NodeId lhs, NodeId rhs) { return lhs.value < rhs.value; }

size_t std::hash<NodeId>::operator()(NodeId id) const noexcept { return id.value; }

bool operator==(EdgeId lhs, EdgeId rhs) { return lhs.value == rhs.value; }

bool operator<(EdgeId lhs, EdgeId rhs) { return lhs.value < rhs.value; }

bool operator==(CladeIdx lhs, CladeIdx rhs) { return lhs.value == rhs.value; }

bool operator<(CladeIdx lhs, CladeIdx rhs) { return lhs.value < rhs.value; }

template <typename Feature, typename View>
auto& GetFeatureStorage(const FeatureReader<Feature, View>* reader) {
  return static_cast<const View&>(*reader).template GetFeatureStorage<Feature>();
}

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
const auto& GetFeature(FeatureReader<Feature, View>* reader) {
  return static_cast<View&>(*reader).template GetFeature<Feature>();
}

template <typename Feature, typename View>
void SetFeature(FeatureWriter<Feature, View>* writer, Feature&& feature) {
  static_cast<View&>(*writer).template SetFeature<Feature>(
      std::forward<Feature>(feature));
}

bool operator==(MutationPosition lhs, MutationPosition rhs) {
  return lhs.value == rhs.value;
}

bool operator<(MutationPosition lhs, MutationPosition rhs) {
  return lhs.value < rhs.value;
}

template <typename T>
EdgeMutations::EdgeMutations(T&& view) : mutations_(view.begin(), view.end()) {}

template <typename View>
const EdgeMutations& FeatureReader<EdgeMutations, View>::GetEdgeMutations() const {
  return GetFeatureStorage(this);
}

template <typename View>
EdgeMutations& FeatureWriter<EdgeMutations, View>::GetEdgeMutations() {
  return GetFeatureStorage(this);
}

template <typename View>
void FeatureWriter<EdgeMutations, View>::SetEdgeMutations(
    EdgeMutations&& edge_mutations) {
  GetFeatureStorage(this) = std::forward<EdgeMutations>(edge_mutations);
}
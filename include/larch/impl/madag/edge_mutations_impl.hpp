bool operator==(MutationPosition lhs, MutationPosition rhs) {
  return lhs.value == rhs.value;
}

bool operator<(MutationPosition lhs, MutationPosition rhs) {
  return lhs.value < rhs.value;
}

template <typename T>
EdgeMutations::EdgeMutations(T&& view) : mutations_{view.begin(), view.end()} {}

EdgeMutations::EdgeMutations(
    std::map<MutationPosition, std::pair<char, char>>&& mutations)
    : mutations_{std::forward<decltype(mutations_)>(mutations)} {}

EdgeMutations EdgeMutations::Copy() const { return EdgeMutations{mutations_}; }

auto EdgeMutations::begin() const -> decltype(mutations_.begin()) {
  return mutations_.begin();
}

auto EdgeMutations::end() const -> decltype(mutations_.end()) {
  return mutations_.end();
}

size_t EdgeMutations::size() const { return mutations_.size(); }

auto EdgeMutations::operator[](MutationPosition pos) -> decltype(mutations_[pos]) {
  Assert(pos.value != NoId);
  return mutations_[pos];
}

auto EdgeMutations::insert(std::pair<MutationPosition, std::pair<char, char>> mut)
    -> decltype(mutations_.insert(mut)) {
  Assert(mut.first.value != NoId);
  return mutations_.insert(mut);
}

bool EdgeMutations::operator==(const EdgeMutations& rhs) const {
  return mutations_ == rhs.mutations_;
}

bool EdgeMutations::operator!=(const EdgeMutations& rhs) const {
  return mutations_ != rhs.mutations_;
}

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

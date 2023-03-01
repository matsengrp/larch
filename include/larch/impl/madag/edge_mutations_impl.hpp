bool operator==(MutationPosition lhs, MutationPosition rhs) {
  return lhs.value == rhs.value;
}

bool operator!=(MutationPosition lhs, MutationPosition rhs) {
  return lhs.value != rhs.value;
}

bool operator<(MutationPosition lhs, MutationPosition rhs) {
  return lhs.value < rhs.value;
}

template <typename T>
EdgeMutations::EdgeMutations(T&& mutations_view,
                             std::enable_if_t<not std::is_same_v<T, EdgeMutations>>*) {
  for (auto mut : mutations_view) {
    mutations_.insert(mut);
  }
}

EdgeMutations::EdgeMutations(
    ContiguousMap<MutationPosition, std::pair<char, char>>&& mutations)
    : mutations_{std::forward<decltype(mutations_)>(mutations)} {}

EdgeMutations EdgeMutations::Copy() const { return EdgeMutations{mutations_.Copy()}; }

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

auto EdgeMutations::insert(std::pair<MutationPosition, std::pair<char, char>> mut) {
  Assert(mut.first.value != NoId);
  return mutations_.insert(mut);
}

bool EdgeMutations::operator==(const EdgeMutations& rhs) const {
  return mutations_ == rhs.mutations_;
}

bool EdgeMutations::operator!=(const EdgeMutations& rhs) const {
  return mutations_ != rhs.mutations_;
}

template <typename CRTP, typename Tag>
const EdgeMutations& FeatureConstView<EdgeMutations, CRTP, Tag>::GetEdgeMutations()
    const {
  return GetFeatureStorage(this);
}

template <typename CRTP, typename Tag>
void FeatureMutableView<EdgeMutations, CRTP, Tag>::SetEdgeMutations(
    EdgeMutations&& edge_mutations) const {
  GetFeatureStorage(this) = std::forward<EdgeMutations>(edge_mutations);
}

template <typename CRTP, typename Tag>
EdgeMutations& FeatureMutableView<EdgeMutations, CRTP, Tag>::GetMutableEdgeMutations()
    const {
  return GetFeatureStorage(this);
}

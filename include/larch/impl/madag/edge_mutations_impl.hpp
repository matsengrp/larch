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
    ContiguousMap<MutationPosition, std::pair<MutationBase, MutationBase>>&& mutations)
    : mutations_{std::forward<decltype(mutations_)>(mutations)} {}

template <typename CRTP>
EdgeMutations EdgeMutations::Copy(const CRTP*) const {
  return EdgeMutations{mutations_.Copy()};
}

auto EdgeMutations::begin() const -> decltype(mutations_.begin()) {
  return mutations_.begin();
}

auto EdgeMutations::end() const -> decltype(mutations_.end()) {
  return mutations_.end();
}

size_t EdgeMutations::size() const { return mutations_.size(); }

bool EdgeMutations::empty() const { return mutations_.empty(); }

auto EdgeMutations::operator[](MutationPosition pos) -> decltype(mutations_[pos]) {
  Assert(pos.value != NoId);
  return mutations_[pos];
}

auto EdgeMutations::insert(
    std::pair<MutationPosition, std::pair<MutationBase, MutationBase>> mut) {
  Assert(mut.first.value != NoId);
  return mutations_.insert(mut);
}

bool EdgeMutations::HasMutationAtPosition(MutationPosition pos) const {
  return mutations_.find(pos) != mutations_.end();
}

std::pair<MutationBase, MutationBase> EdgeMutations::GetMutation(
    MutationPosition pos) const {
  auto it = mutations_.find(pos);
  Assert(it != mutations_.end() and it->first == pos);
  return it->second;
}

bool EdgeMutations::operator==(const EdgeMutations& rhs) const {
  return mutations_ == rhs.mutations_;
}

bool EdgeMutations::operator!=(const EdgeMutations& rhs) const {
  return mutations_ != rhs.mutations_;
}

std::string EdgeMutations::ToString() const {
  std::string result = "<";
  for (auto [pos, muts] : mutations_) {
    result += muts.first;
    result += std::to_string(pos.value);
    result += muts.second;
    result += ",";
  }
  result += ">";
  return result;
}

template <typename CRTP, typename Tag>
const EdgeMutations& FeatureConstView<EdgeMutations, CRTP, Tag>::GetEdgeMutations()
    const {
  return GetFeatureStorage(this).get();
}

template <typename CRTP, typename Tag>
void FeatureMutableView<EdgeMutations, CRTP, Tag>::SetEdgeMutations(
    EdgeMutations&& edge_mutations) const {
  GetFeatureStorage(this).get() = std::forward<EdgeMutations>(edge_mutations);
}

template <typename CRTP, typename Tag>
EdgeMutations& FeatureMutableView<EdgeMutations, CRTP, Tag>::GetMutableEdgeMutations()
    const {
  return GetFeatureStorage(this).get();
}

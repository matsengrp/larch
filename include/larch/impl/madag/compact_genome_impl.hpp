const CompactGenome* CompactGenome::GetEmpty() {
  static const CompactGenome empty = {};
  return &empty;
}

static inline void AssertMut(MutationPosition, [[maybe_unused]] MutationBase mut) {
  Assert(mut == 'A' or mut == 'C' or mut == 'G' or mut == 'T' or mut == 'N');
}

static void ComputeMutations(const EdgeMutations& edge_mutations,
                             std::string_view reference_sequence,
                             ContiguousMap<MutationPosition, MutationBase>& result) {
  for (auto [pos, nucs] : edge_mutations) {
    const bool is_valid = nucs.second != reference_sequence.at(pos.value - 1);
    auto it = result.find(pos);
    if (it != result.end() and it->first == pos) {
      if (is_valid) {
        AssertMut(pos, nucs.second);
        it->second = nucs.second;
      } else {
        result.erase(it);
      }
    } else {
      if (is_valid) {
        AssertMut(pos, nucs.second);
        result.insert_or_assign(pos, nucs.second);
      }
    }
  }
}

CompactGenome::CompactGenome(ContiguousMap<MutationPosition, MutationBase>&& mutations)
    : mutations_{std::move(mutations)}, hash_{ComputeHash(mutations_)} {
#ifndef NDEBUG
  for (auto [pos, mut] : mutations_) {
    AssertMut(pos, mut);
  }
#endif
}

CompactGenome::CompactGenome(ContiguousMap<MutationPosition, MutationBase>&& mutations,
                             size_t hash)
    : mutations_{std::move(mutations)}, hash_{hash} {
#ifndef NDEBUG
  for (auto [pos, mut] : mutations_) {
    AssertMut(pos, mut);
  }
#endif
}

CompactGenome::CompactGenome(const std::string& sequence,
                             const std::string& reference_sequence) {
  Assert(sequence.size() == reference_sequence.size());
  for (size_t i = 0; i < sequence.size(); i++) {
    if (sequence[i] != reference_sequence[i]) {
      AssertMut({i + 1}, {sequence[i]});
      mutations_.insert({{i + 1}, {sequence[i]}});
    }
  }
  hash_ = ComputeHash(mutations_);
}

void CompactGenome::AddParentEdge(const EdgeMutations& mutations,
                                  const CompactGenome& parent,
                                  std::string_view reference_sequence) {
  mutations_.Union(parent.mutations_);
  ComputeMutations(mutations, reference_sequence, mutations_);
  hash_ = ComputeHash(mutations_);
}

void CompactGenome::ApplyChanges(
    const ContiguousMap<MutationPosition, MutationBase>& changes) {
  for (auto change : changes) {
    AssertMut(change.first, change.second);
    mutations_.insert_or_assign(change.first, change.second);
  }
}

bool CompactGenome::HasMutationAtPosition(MutationPosition pos) const {
  auto it = mutations_.find(pos);
  return (it != mutations_.end() and it->first == pos);
}

MutationBase CompactGenome::GetBase(MutationPosition pos,
                                    std::string_view reference_sequence) const {
  auto it = mutations_.find(pos);
  if (it != mutations_.end() and it->first == pos) {
    return it->second;
  }
  return reference_sequence.at(pos.value - 1);
}

ContiguousSet<MutationPosition> CompactGenome::DifferingSites(
    const CompactGenome& other) const {
  ContiguousSet<MutationPosition> result;
  for (auto [pos, base] : mutations_) {
    auto it = other.mutations_.find(pos);
    if (it == other.mutations_.end() or !it->second.IsCompatible(base)) {
      result.insert(pos);
    }
  }
  for (auto [pos, base] : other.mutations_) {
    auto it = mutations_.find(pos);
    if (it == mutations_.end() or !it->second.IsCompatible(base)) {
      result.insert(pos);
    }
  }
  return result;
}

bool CompactGenome::operator==(const CompactGenome& rhs) const noexcept {
  if (hash_ != rhs.hash_) {
    return false;
  }
  return mutations_ == rhs.mutations_;
}

bool CompactGenome::operator!=(const CompactGenome& rhs) const noexcept {
  return mutations_ != rhs.mutations_;
}

bool CompactGenome::operator<(const CompactGenome& rhs) const noexcept {
  return mutations_ < rhs.mutations_;
}

bool CompactGenome::IsCompatible(const CompactGenome& rhs,
                                 std::string_view reference_sequence) const {
  const auto& lhs = *this;
  for (auto [pos, mut] : lhs) {
    if (!rhs.GetBase(pos, reference_sequence)
             .IsCompatible(lhs.GetBase(pos, reference_sequence))) {
      return false;
    }
  }
  for (auto [pos, mut] : rhs) {  // NOLINT(readability-use-anyofallof)
    if (lhs.HasMutationAtPosition(pos)) {
      continue;
    }
    if (!lhs.GetBase(pos, reference_sequence)
             .IsCompatible(rhs.GetBase(pos, reference_sequence))) {
      return false;
    }
  }
  return true;
}

bool CompactGenome::ContainsAmbiguity() const {
  for (auto [pos, base] : mutations_) {  // NOLINT(readability-use-anyofallof)
    if (base.IsAmbiguous()) {
      return true;
    }
  }
  return false;
}

size_t CompactGenome::Hash() const noexcept { return hash_; }

std::optional<MutationBase> CompactGenome::operator[](MutationPosition pos) const {
  auto it = mutations_.find(pos);
  if (it != mutations_.end() and it->first == pos) {
    return it->second;
  }
  return std::nullopt;
}

auto CompactGenome::begin() const -> decltype(mutations_.begin()) {
  return mutations_.begin();
}

auto CompactGenome::end() const -> decltype(mutations_.end()) {
  return mutations_.end();
}

bool CompactGenome::empty() const { return mutations_.empty(); }

template <typename CRTP>
CompactGenome CompactGenome::Copy(const CRTP*) const {
  CompactGenome result{mutations_.Copy(), hash_};
  return result;
}

std::string CompactGenome::ToString() const {
  std::string result = "<";
  for (auto mutpair : mutations_) {
    result += std::to_string(mutpair.first.value);
    result += mutpair.second;
    result += ",";
  }
  result += ">";
  return result;
}

std::string CompactGenome::ToSequence(std::string_view reference_sequence) const {
  std::string sequence(reference_sequence);
  for (auto [pos, base] : mutations_) {
    sequence[pos.value - 1] = base.ToChar();
  }
  return sequence;
}

EdgeMutations CompactGenome::ToEdgeMutations(std::string_view reference_sequence,
                                             const CompactGenome& parent,
                                             const CompactGenome& child) {
  EdgeMutations result;
  for (auto [pos, child_base] : child) {
    MutationBase parent_base = reference_sequence.at(pos.value - 1);
    auto opt_parent_base = parent[pos];
    if (opt_parent_base.has_value()) {
      parent_base = opt_parent_base.value();
    }
    if (!parent_base.IsCompatible(child_base)) {
      result[pos] = {parent_base, child_base.GetFirstBase()};
    }
  }

  for (auto [pos, parent_base] : parent) {
    MutationBase child_base = reference_sequence.at(pos.value - 1);
    auto opt_child_base = child[pos];
    if (opt_child_base.has_value()) {
      child_base = opt_child_base.value();
    }
    if (!child_base.IsCompatible(parent_base)) {
      result[pos] = {parent_base, child_base.GetFirstBase()};
    }
  }
  return result;
}

size_t CompactGenome::ComputeHash(
    const ContiguousMap<MutationPosition, MutationBase>& mutations) {
  size_t result = 0;
  for (auto [pos, base] : mutations) {
    result = HashCombine(result, pos.value);
    result = HashCombine(result, static_cast<size_t>(base.ToChar()));
  }
  return result;
}

std::size_t std::hash<CompactGenome>::operator()(
    const CompactGenome& cg) const noexcept {
  return cg.Hash();
}

bool std::equal_to<CompactGenome>::operator()(const CompactGenome& lhs,
                                              const CompactGenome& rhs) const noexcept {
  return lhs == rhs;
}

template <typename CRTP, typename Tag>
const CompactGenome& FeatureConstView<CompactGenome, CRTP, Tag>::GetCompactGenome()
    const {
  return GetFeatureStorage(this).get();
}

template <typename CRTP, typename Tag>
// NOLINTNEXTLINE(cppcoreguidelines-c-copy-assignment-signature,misc-unconventional-assign-operator)
auto& FeatureMutableView<CompactGenome, CRTP, Tag>::operator=(
    CompactGenome&& compact_genome) const {
  GetFeatureStorage(this).get() = std::forward<CompactGenome>(compact_genome);
  return *this;
}

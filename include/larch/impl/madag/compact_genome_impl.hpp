const CompactGenome* CompactGenome::Empty() {
  static const CompactGenome empty = {};
  return &empty;
}

static inline void AssertMut(MutationPosition, MutationBase mut) {
  Assert(mut == 'A' or mut == 'C' or mut == 'G' or mut == 'T');
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
        result.Insert(it, {pos, nucs.second});
      }
    }
  }
}

CompactGenome::CompactGenome(ContiguousMap<MutationPosition, MutationBase>&& mutations)
    : mutations_{std::move(mutations)}, hash_{ComputeHash(mutations_)} {
  for (auto [pos, mut] : mutations_) {
    AssertMut(pos, mut);
  }
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
    mutations_.insert(change);
  }
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
    if (it == other.mutations_.end() or it->second != base) {
      result.insert(pos);
    }
  }
  for (auto [pos, base] : other.mutations_) {
    auto it = mutations_.find(pos);
    if (it == mutations_.end() or it->second != base) {
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

CompactGenome CompactGenome::Copy() const {
  CompactGenome result{mutations_.Copy()};
  result.hash_ = hash_;
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
    if (parent_base != child_base) {
      result[pos] = {parent_base, child_base};
    }
  }

  for (auto [pos, parent_base] : parent) {
    MutationBase child_base = reference_sequence.at(pos.value - 1);
    auto opt_child_base = child[pos];
    if (opt_child_base.has_value()) {
      child_base = opt_child_base.value();
    }
    if (child_base != parent_base) {
      result[pos] = {parent_base, child_base};
    }
  }
  return result;
}

size_t CompactGenome::ComputeHash(
    const ContiguousMap<MutationPosition, MutationBase>& mutations) {
  size_t result = 0;
  for (auto [pos, base] : mutations) {
    result = HashCombine(result, pos.value);
    result = HashCombine(result, base.ToChar());
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
  return GetFeatureStorage(this);
}

template <typename CRTP, typename Tag>
auto& FeatureMutableView<CompactGenome, CRTP, Tag>::operator=(
    CompactGenome&& compact_genome) const {
  GetFeatureStorage(this) = std::forward<CompactGenome>(compact_genome);
  return *this;
}

const CompactGenome* CompactGenome::Empty() {
  static CompactGenome empty = {};
  return &empty;
}

static void ComputeMutations(const EdgeMutations& edge_mutations,
                             std::string_view reference_sequence,
                             std::vector<std::pair<MutationPosition, char>>& result) {
  for (auto [pos, nucs] : edge_mutations) {
    const bool is_valid = nucs.second != reference_sequence.at(pos.value - 1);
    auto it = std::lower_bound(result.begin(), result.end(), pos,
                               [](std::pair<MutationPosition, char> lhs,
                                  MutationPosition rhs) { return lhs.first < rhs; });
    if (it != result.end() and it->first == pos) {
      if (is_valid) {
        it->second = nucs.second;
      } else {
        result.erase(it);
      }
    } else {
      if (is_valid) {
        result.insert(it, {pos, nucs.second});
      }
    }
  }
}

CompactGenome::CompactGenome(std::vector<std::pair<MutationPosition, char>>&& mutations)
    : mutations_{mutations}, hash_{ComputeHash(mutations_)} {}

void CompactGenome::AddParentEdge(const EdgeMutations& mutations,
                                  const CompactGenome& parent,
                                  std::string_view reference_sequence) {
  std::vector<std::pair<MutationPosition, char>> mutations_copy = mutations_;
  mutations_.clear();
  std::set_union(mutations_copy.begin(), mutations_copy.end(),
                 parent.mutations_.begin(), parent.mutations_.end(),
                 std::back_inserter(mutations_));
  ComputeMutations(mutations, reference_sequence, mutations_);
  hash_ = ComputeHash(mutations_);
}

bool CompactGenome::operator==(const CompactGenome& rhs) const noexcept {
  if (hash_ != rhs.hash_) {
    return false;
  }
  return mutations_ == rhs.mutations_;
}

bool CompactGenome::operator<(const CompactGenome& rhs) const noexcept {
  return mutations_ < rhs.mutations_;
}

size_t CompactGenome::Hash() const noexcept { return hash_; }

std::optional<char> CompactGenome::operator[](MutationPosition pos) const {
  auto it = std::lower_bound(mutations_.begin(), mutations_.end(), pos,
                             [](std::pair<MutationPosition, char> lhs,
                                MutationPosition rhs) { return lhs.first < rhs; });
  if (it != mutations_.end() and it->first == pos) {
    return it->second;
  } else {
    return std::nullopt;
  }
}

auto CompactGenome::begin() const -> decltype(mutations_.begin()) {
  return mutations_.begin();
}

auto CompactGenome::end() const -> decltype(mutations_.end()) {
  return mutations_.end();
}

bool CompactGenome::empty() const { return mutations_.empty(); }

CompactGenome CompactGenome::Copy() const {
  CompactGenome result;
  result.mutations_ = mutations_;
  result.hash_ = hash_;
  return result;
}

std::string CompactGenome::ToString() const {
  std::string result = "\n<";
  for (auto mutpair : mutations_) {
    result += std::to_string(mutpair.first.value);
    result += mutpair.second;
    result += ",";
  }
  result += ">\n";
  return result;
}

EdgeMutations CompactGenome::ToEdgeMutations(std::string_view reference_sequence,
                                             const CompactGenome& parent,
                                             const CompactGenome& child) {
  EdgeMutations result;
  for (auto [pos, child_base] : child) {
    char parent_base = reference_sequence.at(pos.value - 1);
    auto opt_parent_base = parent[pos];
    if (opt_parent_base.has_value()) {
      parent_base = opt_parent_base.value();
    }
    if (parent_base != child_base) {
      result[pos] = {parent_base, child_base};
    }
  }

  for (auto [pos, parent_base] : parent) {
    char child_base = reference_sequence.at(pos.value - 1);
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
    const std::vector<std::pair<MutationPosition, char>>& mutations) {
  size_t result = 0;
  for (auto [pos, base] : mutations) {
    result = HashCombine(result, pos.value);
    result = HashCombine(result, static_cast<size_t>(base));
  }
  return result;
}

template <typename View>
const CompactGenome& FeatureReader<CompactGenome, View>::GetCompactGenome() const {
  return GetFeatureStorage(this);
}

template <typename View>
void FeatureWriter<CompactGenome, View>::SetCompactGenome(
    CompactGenome&& compact_genome) {
  GetFeatureStorage(this) = std::forward<CompactGenome>(compact_genome);
}

std::size_t std::hash<CompactGenome>::operator()(
    const CompactGenome& cg) const noexcept {
  return cg.Hash();
}

std::size_t std::equal_to<CompactGenome>::operator()(
    const CompactGenome& lhs, const CompactGenome& rhs) const noexcept {
  return lhs == rhs;
}

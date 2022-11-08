// CladeUnion::CladeUnion(const std::vector<const CladeUnion*>& clades)
//     : leafs_{clades |
//              ranges::views::transform([](const CladeUnion* i) { return i.leafs_; }) |
//              ranges::actions::join | ranges::actions::sort},
//       hash_{ComputeHash(leafs_)} {}

CladeUnion::CladeUnion(const std::vector<const CompactGenome*>& leafs)
    : leafs_{leafs}, hash_{ComputeHash(leafs_)} {}

bool CladeUnion::operator==(const CladeUnion& rhs) const noexcept {
  return leafs_ == rhs.leafs_;
}

size_t CladeUnion::Hash() const noexcept { return hash_; }

auto CladeUnion::begin() const { return leafs_.begin(); }

auto CladeUnion::end() const { return leafs_.end(); }

bool CladeUnion::empty() const { return leafs_.empty(); }

size_t CladeUnion::size() const { return leafs_.size(); }

size_t CladeUnion::ComputeHash(
    const std::vector<const CompactGenome*>& leafs) noexcept {
  size_t hash = 0;
  for (const auto* leaf : leafs) {
    hash = HashCombine(hash, leaf->Hash());
  }
  return hash;
}

std::size_t std::hash<CladeUnion>::operator()(const CladeUnion& cu) const noexcept {
  return cu.Hash();
}

bool std::equal_to<CladeUnion>::operator()(const CladeUnion& lhs,
                                           const CladeUnion& rhs) const noexcept {
  return lhs == rhs;
}
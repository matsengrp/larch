
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

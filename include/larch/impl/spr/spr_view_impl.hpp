bool FitchSet::operator==(const FitchSet&) const noexcept { return true; }

bool FitchSet::operator<(const FitchSet&) const noexcept { return false; }

size_t FitchSet::Hash() const noexcept { return 0; }

FitchSet FitchSet::Copy() const {
  FitchSet result;
  return result;
}

std::size_t std::hash<FitchSet>::operator()(const FitchSet& fs) const noexcept {
  return fs.Hash();
}

bool std::equal_to<FitchSet>::operator()(const FitchSet& lhs,
                                         const FitchSet& rhs) const noexcept {
  return lhs == rhs;
}

template <typename CRTP, typename Tag>
const FitchSet& FeatureConstView<FitchSet, CRTP, Tag>::GetFitchSet() const {
  return GetFeatureStorage(this);
}

template <typename CRTP, typename Tag>
auto& FeatureMutableView<FitchSet, CRTP, Tag>::operator=(FitchSet&& fitch_set) {
  GetFeatureStorage(this) = std::forward<FitchSet>(fitch_set);
  return *this;
}
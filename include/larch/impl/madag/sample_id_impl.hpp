
const SampleId* SampleId::GetEmpty() {
  static const SampleId empty = {};
  return &empty;
}

std::size_t std::hash<SampleId>::operator()(const SampleId& sid) const noexcept {
  return sid.Hash();
}

bool std::equal_to<SampleId>::operator()(const SampleId& lhs,
                                         const SampleId& rhs) const noexcept {
  if (lhs.hash_ != rhs.hash_) {
    return false;
  }
  return lhs.sample_id_ == rhs.sample_id_;
}

template <typename CRTP, typename Tag>
const std::optional<std::string>& FeatureConstView<SampleId, CRTP, Tag>::GetSampleId()
    const {
  return GetFeatureStorage(this).get().sample_id_;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<SampleId, CRTP, Tag>::SetSampleId(
    const std::optional<std::string>& sample_id) const {
  GetFeatureStorage(this).get().hash_ = SampleId::ComputeHash(sample_id);
  GetFeatureStorage(this).get().sample_id_ = sample_id;
}

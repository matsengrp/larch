
std::size_t std::hash<SampleId>::operator()(const SampleId &sid) const noexcept {
  if (not sid.sample_id_.has_value()) {
    return 0;
  }
  return std::hash<std::string>{}(sid.sample_id_.value());
}

bool std::equal_to<SampleId>::operator()(const SampleId &lhs,
                                         const SampleId &rhs) const noexcept {
  return lhs.sample_id_ == rhs.sample_id_;
}

template <typename CRTP, typename Tag>
const std::optional<std::string> &FeatureConstView<SampleId, CRTP, Tag>::GetSampleId()
    const {
  return GetFeatureStorage(this).sample_id_;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<SampleId, CRTP, Tag>::SetSampleId(
    const std::optional<std::string> &sample_id) const {
  GetFeatureStorage(this).sample_id_ = sample_id;
}

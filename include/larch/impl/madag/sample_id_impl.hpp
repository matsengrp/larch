std::size_t std::hash<SampleIdStorage>::operator()(
    const SampleIdStorage& sid) const noexcept {
  return sid.Hash();
}

bool std::equal_to<SampleIdStorage>::operator()(
    const SampleIdStorage& lhs, const SampleIdStorage& rhs) const noexcept {
  return lhs.Value() == rhs.Value();
}

const SampleId* SampleId::GetEmpty() {
  static const SampleId empty = {};
  return &empty;
}

std::size_t std::hash<SampleId>::operator()(const SampleId& sid) const noexcept {
  return sid.Hash();
}

bool std::equal_to<SampleId>::operator()(const SampleId& lhs,
                                         const SampleId& rhs) const noexcept {
  return lhs.target_ == rhs.target_;
}

template <typename CRTP, typename Tag>
std::optional<std::string_view> FeatureConstView<SampleId, CRTP, Tag>::GetSampleId()
    const {
  auto* target_ = GetFeatureStorage(this).get().target_;
  if (not target_) {
    return std::nullopt;
  } else {
    return target_->Value();
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<SampleId, CRTP, Tag>::SetSampleId(
    std::optional<std::string_view> sample_id) const {
  if (not sample_id.has_value()) {
    GetFeatureStorage(this).get().target_ = nullptr;
  } else {
    GetFeatureStorage(this).get().target_ =
        &SampleIdStorage::Get(std::string{sample_id.value()});
  }
}

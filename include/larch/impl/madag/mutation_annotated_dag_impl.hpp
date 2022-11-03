template <typename View>
const std::string& FeatureReader<ReferenceSequence, View>::GetReferenceSequence()
    const {
  return GetFeatureStorage(this).reference_sequence_;
}

template <typename View>
void FeatureWriter<ReferenceSequence, View>::SetReferenceSequence(
    std::string_view reference_sequence) const {
  GetFeatureStorage(this).reference_sequence_ = reference_sequence;
}
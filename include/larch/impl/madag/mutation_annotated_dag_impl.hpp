template <typename View>
const std::string& FeatureReader<ReferenceSequence, View>::GetReferenceSequence()
    const {
  return GetFeatureStorage(this).reference_sequence_;
}

template <typename View>
void FeatureReader<ReferenceSequence, View>::AssertUA() const {
  const View& dag = static_cast<const View&>(*this);
  Assert(dag.HaveRoot());
  typename View::Node ua = dag.GetRoot();
  Assert(ua.GetCladesCount() == 1);
}

template <typename View>
bool FeatureReader<ReferenceSequence, View>::HaveUA() const {
  const View& dag = static_cast<const View&>(*this);
  Assert(dag.HaveRoot());
  typename View::Node ua = dag.GetRoot();
  if (ua.GetCladesCount() != 1) {
    return false;
  }
  return true;
}

template <typename View>
void FeatureWriter<ReferenceSequence, View>::SetReferenceSequence(
    std::string_view reference_sequence) const {
  GetFeatureStorage(this).reference_sequence_ = reference_sequence;
}

template <typename View>
void FeatureWriter<ReferenceSequence, View>::AddUA(
    const EdgeMutations& mutations_at_root) const {
  const View& dag = static_cast<const View&>(*this);
  Assert(not dag.HaveUA());
  typename View::Node root = dag.GetRoot();
  typename View::Node ua_node = dag.AppendNode();
  typename View::Edge ua_edge = dag.AppendEdge(ua_node, root, {0});
  ua_edge.SetEdgeMutations(mutations_at_root.Copy());
  dag.BuildConnections();
  dag.AssertUA();
}

template <typename View>
const std::optional<std::string>& FeatureReader<SampleId, View>::GetSampleId() const {
  return GetFeatureStorage(this).sample_id_;
}

template <typename View>
void FeatureWriter<SampleId, View>::SetSampleId(
    const std::optional<std::string>& sample_id) const {
  GetFeatureStorage(this).sample_id_ = sample_id;
}
#include <vector>
template <typename View>
const std::string& FeatureReader<ReferenceSequence, View>::GetReferenceSequence() {
  return GetFeature(this).reference_sequence_;
}

template <typename View>
void FeatureReader<ReferenceSequence, View>::AssertUA() {
  View& dag = static_cast<View&>(*this);
  Assert(dag.HaveRoot());
  typename View::Node ua = dag.GetRoot();
  Assert(ua.GetCladesCount() == 1);
}

template <typename View>
bool FeatureReader<ReferenceSequence, View>::HaveUA() {
  View& dag = static_cast<View&>(*this);
  Assert(dag.HaveRoot());
  typename View::Node ua = dag.GetRoot();
  return ua.GetCladesCount() == 1;
}

template <typename View>
void FeatureWriter<ReferenceSequence, View>::SetReferenceSequence(
    std::string_view reference_sequence) {
  GetFeature(this).reference_sequence_ = reference_sequence;
}

template <typename View>
void FeatureWriter<ReferenceSequence, View>::AddUA(
    const EdgeMutations& mutations_at_root) {
  View& dag = static_cast<View&>(*this);
  Assert(not dag.HaveUA());
  typename View::Node root = dag.GetRoot();
  typename View::Node ua_node = dag.AppendNode();
  typename View::Edge ua_edge = dag.AppendEdge(ua_node, root, {0});
  ua_edge.SetEdgeMutations(mutations_at_root.Copy());
  dag.BuildConnections();
  dag.AssertUA();
}

template <typename View>
void FeatureWriter<ReferenceSequence, View>::RecomputeCompactGenomes() {
  View& dag = static_cast<View&>(*this);
  using Node = typename View::Node;
  using Edge = typename View::Edge;

  std::vector<CompactGenome> result_cgs;
  result_cgs.resize(dag.GetNodesCount());

  auto ComputeCG = [&](auto& self, Node for_node, std::vector<CompactGenome>& result) {
    CompactGenome& compact_genome = result.at(for_node.GetId().value);
    if (for_node.IsRoot()) {
      compact_genome = {};
      return;
    }
    if (not compact_genome.empty()) {
      return;
    }
    Edge edge = *(for_node.GetParents().begin());
    self(self, edge.GetParent(), result);
    const EdgeMutations& mutations = edge.GetEdgeMutations();
    const CompactGenome& parent = result.at(edge.GetParentId().value);
    compact_genome.AddParentEdge(mutations, parent, dag.GetReferenceSequence());
  };
  for (Node node : dag.GetNodes()) {
    ComputeCG(ComputeCG, node, result_cgs);
  }
  for (Node node : dag.GetNodes()) {
    node.SetCompactGenome(std::move(result_cgs.at(node.GetId().value)));
  }
  std::unordered_map<CompactGenome, NodeId> leaf_cgs;
  for (Node node : dag.GetNodes()) {
    if (node.IsLeaf()) {
      bool success =
          leaf_cgs.emplace(node.GetCompactGenome().Copy(), node.GetId()).second;
      if (not success) {
        // std::cout << "Error in ComputeCompactGenomes: had a non-unique leaf node at "
        //           << node.GetId().value << " also seen at "
        //           << leaf_cgs.at(node.GetCompactGenome().Copy()).value
        //           << "\nCompact Genome is\n"
        //           << node.GetCompactGenome().ToString() << "\n"
        //           << std::flush;
        Fail("Error in ComputeCompactGenomes: had a non-unique leaf node");
      }
    }
  }
}

template <typename View>
void FeatureWriter<ReferenceSequence, View>::RecomputeEdgeMutations() {
  View& dag = static_cast<View&>(*this);
  using Edge = typename View::Edge;

  for (Edge edge : dag.GetEdges()) {
    edge.SetEdgeMutations(CompactGenome::ToEdgeMutations(
        dag.GetReferenceSequence(), edge.GetParent().GetCompactGenome(),
        edge.GetChild().GetCompactGenome()));
  }
}

template <typename View>
const std::optional<std::string>& FeatureReader<SampleId, View>::GetSampleId() {
  return GetFeature(this).sample_id_;
}

template <typename View>
void FeatureWriter<SampleId, View>::SetSampleId(
    std::optional<std::string>&& sample_id) {
  GetFeature(this).sample_id_ = sample_id;
}

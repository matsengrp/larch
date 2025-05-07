template <typename CRTP, typename Tag>
const std::string&
FeatureConstView<ReferenceSequence, CRTP, Tag>::GetReferenceSequence() const {
  return GetFeatureStorage(this).get().reference_sequence_;
}

template <typename CRTP, typename Tag>
void FeatureConstView<ReferenceSequence, CRTP, Tag>::AssertUA() const {
  [[maybe_unused]] auto& dag = static_cast<const CRTP&>(*this);
  Assert(dag.HaveRoot());
  // auto ua = dag.GetRoot();
  // TODO Assert((not dag.IsTree()) or ua.GetCladesCount() == 1);
}

template <typename CRTP, typename Tag>
bool FeatureConstView<ReferenceSequence, CRTP, Tag>::HaveUA() const {
  auto& dag = static_cast<const CRTP&>(*this);
  Assert(dag.HaveRoot());
  auto ua = dag.GetRoot();
  return ua.GetCladesCount() == 1;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<ReferenceSequence, CRTP, Tag>::SetReferenceSequence(
    std::string_view reference_sequence) const {
  GetFeatureStorage(this).get().reference_sequence_ = reference_sequence;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<ReferenceSequence, CRTP, Tag>::
    SetCompactGenomesFromNodeSequenceMap(const NodeSeqMap& node_sequence_map) const {
  auto dag = static_cast<const CRTP&>(*this);
  using Node = typename decltype(dag)::NodeView;

  auto ref_seq = dag.GetReferenceSequence();
  auto ComputeCGFromSequence = [&ref_seq](const std::string& leaf_seq, Node for_node) {
    CompactGenome new_cg(leaf_seq, ref_seq);
    for_node = std::move(new_cg);
  };

  for (auto node : dag.GetNodes()) {
    if (node_sequence_map.find(node.GetId()) != node_sequence_map.end()) {
      auto& node_seq = node_sequence_map.find(node.GetId())->second;
      ComputeCGFromSequence(node_seq, node);
    }
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<ReferenceSequence, CRTP, Tag>::
    SetCompactGenomesFromNodeMutationMap(NodeMutMap&& node_mutation_map) const {
  auto dag = static_cast<const CRTP&>(*this);
  using Node = typename decltype(dag)::NodeView;

  auto ComputeCGFromMutation =
      [&dag](ContiguousMap<MutationPosition, MutationBase>&& new_muts, Node for_node) {
        CompactGenome new_cg(std::move(new_muts));
        for_node = std::move(new_cg);
      };

  for (auto node : dag.GetNodes()) {
    if (node_mutation_map.find(node.GetId()) != node_mutation_map.end()) {
      auto&& muts = node_mutation_map.find(node.GetId())->second;
      ComputeCGFromMutation(std::move(muts), node);
    }
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<ReferenceSequence, CRTP, Tag>::
    UpdateCompactGenomesFromNodeMutationMap(NodeMutMap&& node_mutation_map) const {
  auto dag = static_cast<const CRTP&>(*this);
  using Node = typename decltype(dag)::NodeView;

  auto ComputeCGFromMutation =
      [this](ContiguousMap<MutationPosition, MutationBase>&& new_muts, Node for_node) {
        CompactGenome new_cg =
            for_node.GetCompactGenome().Copy(static_cast<const CRTP*>(this));
        new_cg.ApplyChanges(new_muts);
        for_node = std::move(new_cg);
      };

  for (auto node : dag.GetNodes()) {
    if (node_mutation_map.find(node.GetId()) != node_mutation_map.end()) {
      auto&& muts = node_mutation_map.find(node.GetId())->second;
      ComputeCGFromMutation(std::move(muts), node);
    }
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<ReferenceSequence, CRTP, Tag>::AddUA(
    const EdgeMutations& mutations_at_root) const {
  auto dag = static_cast<const CRTP&>(*this);
  using Node = typename decltype(dag)::NodeView;
  using Edge = typename decltype(dag)::EdgeView;

  // Assert(not dag.HaveUA());
  Node root = dag.GetRoot();
  Node ua_node = dag.AppendNode();
  Edge ua_edge = dag.AppendEdge(ua_node, root, {0});
  ua_edge.SetEdgeMutations(mutations_at_root.Copy(static_cast<const CRTP*>(this)));
  dag.BuildConnections();
  dag.AssertUA();
}

template <typename CRTP, typename Tag>
template <IdContinuity Cont>
void FeatureMutableView<ReferenceSequence, CRTP, Tag>::RecomputeCompactGenomes(
    bool recompute_leaves) const {
  auto dag = static_cast<const CRTP&>(*this);
  using Node = typename decltype(dag)::NodeView;
  using Edge = typename decltype(dag)::EdgeView;

  auto max_id = dag.GetNodesCount();
  for (auto node : dag.GetNodes()) {
    max_id = max_id > node.GetId().value ? max_id : node.GetId().value;
  }
  IdContainer<NodeId, CompactGenome, Cont> new_cgs;
  new_cgs.resize(max_id + 1);
  auto ComputeCG = [&new_cgs, dag](auto& self, Node for_node) {
    CompactGenome& compact_genome = new_cgs.at(for_node.GetId());
    if (for_node.IsUA()) {
      compact_genome = {};
      return;
    }
    if (not compact_genome.empty()) {
      return;
    }
    Edge edge = for_node.GetFirstParent();
    self(self, edge.GetParent());
    const EdgeMutations& mutations = edge.GetEdgeMutations();
    const CompactGenome& parent = new_cgs.at(edge.GetParentId());
    compact_genome.AddParentEdge(mutations, parent, dag.GetReferenceSequence());
  };
  for (Node node : dag.GetNodes()) {
    if (recompute_leaves || !node.IsLeaf()) {
      ComputeCG(ComputeCG, node);
    }
  }
  for (Node node : dag.GetNodes()) {
    if (recompute_leaves || !node.IsLeaf()) {
      node = std::move(new_cgs.at(node.GetId()));
    }
  }
#ifndef NDEBUG
  // TODO extract validation to separate function to not hurt performance
  ContiguousMap<CompactGenome, NodeId> leaf_cgs;
  for (Node node : dag.GetNodes()) {
    if (node.IsLeaf()) {
      bool success =
          leaf_cgs
              .insert({node.GetCompactGenome().Copy(static_cast<const CRTP*>(this)),
                       node.GetId()})
              .second;
      if (not success) {
        // std::cout << "Error in ComputeCompactGenomes: had a non-unique leaf node at "
        //           << node.GetId().value << " also seen at "
        //           << leaf_cgs.at(node.GetCompactGenome().Copy()).value
        //           << "\nCompact Genome is\n"
        //           << node.GetCompactGenome().ToString() << "\n"
        //           << std::flush;
        // TODO Fail("Error in ComputeCompactGenomes: had a non-unique leaf node");
      }
    }
  }
#endif
}

template <typename CRTP, typename Tag>
void FeatureMutableView<ReferenceSequence, CRTP, Tag>::SampleIdsFromCG(
    bool coerce) const {
  auto dag = static_cast<const CRTP&>(*this);
  for (auto leaf : dag.GetLeafs()) {
    if (not leaf.HaveSampleId() or coerce) {
      std::string id = leaf.GetCompactGenome().ToString();
      Assert(not id.empty());
      if constexpr (decltype(leaf)::template contains_feature<Deduplicate<SampleId>>) {
        auto id_iter = dag.template AsFeature<Deduplicate<SampleId>>().AddDeduplicated(
            SampleId{id});
        leaf = id_iter.first;
      } else {
        leaf.SetSampleId(id);
      }
    }
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<ReferenceSequence, CRTP, Tag>::RecomputeEdgeMutations() const {
  auto dag = static_cast<const CRTP&>(*this);
  using Edge = typename decltype(dag)::EdgeView;

  for (Edge edge : dag.GetEdges()) {
    edge.SetEdgeMutations(CompactGenome::ToEdgeMutations(
        dag.GetReferenceSequence(), edge.GetParent().GetCompactGenome(),
        edge.GetChild().GetCompactGenome()));
  }
}

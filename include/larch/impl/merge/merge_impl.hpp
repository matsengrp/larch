Merge::Merge(std::string_view reference_sequence) {
  result_dag_storage_.View().SetReferenceSequence(reference_sequence);
}

void Merge::AddDAGs(const std::vector<MADAG>& dags) {
  for (auto dag : dags) {
    dag.AssertUA();
  }
  std::vector<size_t> dag_idxs;
  dag_idxs.resize(dags.size());
  std::iota(dag_idxs.begin(), dag_idxs.end(), dags_.size());

  for (MADAG i : dags) {
    dags_.emplace_back(i);
  }
  dags_node_labels_.resize(dags_.size());

  tbb::parallel_for_each(dag_idxs.begin(), dag_idxs.end(), [&](size_t dag_idx) {
    MADAG dag = dags_.at(dag_idx);
    std::vector<NodeLabel>& labels = dags_node_labels_.at(dag_idx);
    labels.resize(dag.GetNodesCount());
    for (size_t node_idx = 0; node_idx < dag.GetNodesCount(); ++node_idx) {
      const CompactGenome* cg =
          result_dag_storage_.View()
              .GetNodesContainer()
              .GetFeatureGlobalData<Deduplicate<CompactGenome>>()
              .AddUnique(dag.Get(NodeId{node_idx}).GetCompactGenome().Copy());
      labels.at(node_idx).SetCompactGenome(cg);
    }
  });

  ComputeLeafSets(dag_idxs);
  MergeDAGs(dag_idxs);
  Assert(result_nodes_.size() == result_dag_storage_.View().GetNodesCount());
  Assert(result_edges_.size() == result_dag_storage_.View().GetEdgesCount());
  result_dag_storage_.View().BuildConnections();
}

MergeDAG Merge::GetResult() const { return result_dag_storage_.View(); }

const std::unordered_map<NodeLabel, NodeId>& Merge::GetResultNodes() const {
  return result_nodes_;
}

void Merge::ComputeResultEdgeMutations() {
  for (auto& [label, edge_id] : result_edges_) {
    Assert(label.GetParent().GetCompactGenome());
    Assert(label.GetChild().GetCompactGenome());
    const CompactGenome& parent = *label.GetParent().GetCompactGenome();
    const CompactGenome& child = *label.GetChild().GetCompactGenome();
    result_dag_storage_.View().Get(edge_id).SetEdgeMutations(
        CompactGenome::ToEdgeMutations(
            result_dag_storage_.View().GetReferenceSequence(), parent, child));
  }
}

bool Merge::ContainsLeafset(const LeafSet& leafset) const {
  return result_dag_storage_.View()
      .GetNodesContainer()
      .GetFeatureGlobalData<Deduplicate<LeafSet>>()
      .ContainsUnique(leafset);
}

void Merge::ComputeLeafSets(const std::vector<size_t>& dag_idxs) {
  tbb::parallel_for_each(dag_idxs.begin(), dag_idxs.end(), [&](size_t dag_idx) {
    MADAG dag = dags_.at(dag_idx);
    std::vector<NodeLabel>& labels = dags_node_labels_.at(dag_idx);
    std::vector<LeafSet> computed_ls = ComputeLeafSets(dag, labels);
    for (Node node : dag.GetNodes()) {
      const LeafSet* ls = result_dag_storage_.View()
                              .GetNodesContainer()
                              .GetFeatureGlobalData<Deduplicate<LeafSet>>()
                              .AddUnique(std::move(computed_ls.at(node.GetId().value)));
      labels.at(node.GetId().value).SetLeafSet(ls);
    }
  });
}

void Merge::MergeDAGs(const std::vector<size_t>& dag_idxs) {
  NodeId node_id{result_dag_storage_.View().GetNodesCount()};
  std::mutex mtx;
  tbb::parallel_for_each(dag_idxs.begin(), dag_idxs.end(), [&](size_t dag_idx) {
    const std::vector<NodeLabel>& labels = dags_node_labels_.at(dag_idx);
    for (auto label : labels) {
      std::unique_lock<std::mutex> lock{mtx};
      if (result_nodes_.try_emplace(label, node_id).second) {
        result_dag_storage_.View().AddNode(node_id);
        result_dag_storage_.View().Get(node_id).SetNodeLabel(NodeLabel{label});
        ++node_id.value;
      }
    }
  });
  tbb::concurrent_vector<EdgeLabel> added_edges;
  tbb::parallel_for_each(dag_idxs.begin(), dag_idxs.end(), [&](size_t dag_idx) {
    MADAG dag = dags_.at(dag_idx);
    const std::vector<NodeLabel>& labels = dags_node_labels_.at(dag_idx);
    for (Edge edge : dag.GetEdges()) {
      const auto& parent_label = labels.at(edge.GetParentId().value);
      const auto& child_label = labels.at(edge.GetChildId().value);
      auto ins = result_edges_.insert({{parent_label, child_label}, {}});
      if (ins.second) {
        added_edges.push_back({parent_label, child_label});
      }
    }
  });
  EdgeId edge_id{result_dag_storage_.View().GetEdgesCount()};
  for (auto edge : added_edges) {
    auto parent = result_nodes_.find(edge.GetParent());
    auto child = result_nodes_.find(edge.GetChild());
    Assert(parent != result_nodes_.end());
    Assert(child != result_nodes_.end());
    Assert(parent->second.value < result_dag_storage_.View().GetNodesCount());
    Assert(child->second.value < result_dag_storage_.View().GetNodesCount());
    result_dag_storage_.View().AddEdge(edge_id, parent->second, child->second,
                                       edge.ComputeCladeIdx());
    auto result_edge_it = result_edges_.find(edge);
    Assert(result_edge_it != result_edges_.end());
    result_edge_it->second = edge_id;
    edge_id.value++;
  }
}

std::vector<LeafSet> Merge::ComputeLeafSets(MADAG dag,
                                            const std::vector<NodeLabel>& labels) {
  std::vector<LeafSet> result;
  result.resize(dag.GetNodesCount());
  auto ComputeLS = [&](auto& self, Node for_node) {
    const LeafSet& leaf_set = result.at(for_node.GetId().value);
    if (not leaf_set.empty()) {
      return;
    }
    for (Node child : for_node.GetChildren() | Transform::GetChild()) {
      self(self, child);
    }
    result.at(for_node.GetId().value) = LeafSet{for_node, labels, result};
  };
  for (Node node : dag.GetNodes()) {
    ComputeLS(ComputeLS, node);
  }
  return result;
}

Merge::Merge(std::string_view reference_sequence) : result_dag_{result_dag_storage_} {
  result_dag_.SetReferenceSequence(reference_sequence);
}

void Merge::AddDAGs(const std::vector<MADAG>& dags) {
  for (auto dag : dags) {
    dag.AssertUA();
  }
  std::vector<size_t> dag_idxs;
  dag_idxs.resize(dags.size());
  std::iota(dag_idxs.begin(), dag_idxs.end(), trees_.size());

  for (MADAG i : dags) {
    trees_.emplace_back(i);
  }
  tree_labels_.resize(trees_.size());

  tbb::parallel_for_each(dag_idxs.begin(), dag_idxs.end(), [&](size_t dag_idx) {
    MADAG dag = trees_.at(dag_idx);
    std::vector<NodeLabel>& labels = tree_labels_.at(dag_idx);
    labels.resize(dag.GetNodesCount());
    for (size_t node_idx = 0; node_idx < dag.GetNodesCount(); ++node_idx) {
      auto cg_iter = all_compact_genomes_.insert(
          dag.Get(NodeId{node_idx}).GetCompactGenome().Copy());
      labels.at(node_idx).SetCompactGenome(std::addressof(*cg_iter.first));
    }
  });

  ComputeLeafSets(dag_idxs);
  MergeTrees(dag_idxs);
  Assert(result_nodes_.size() == result_dag_.GetNodesCount());
  Assert(result_edges_.size() == result_dag_.GetEdgesCount());
  result_dag_.BuildConnections();
}

MergeDAG Merge::GetResult() const { return result_dag_; }

const std::unordered_map<NodeLabel, NodeId>& Merge::GetResultNodes() const {
  return result_nodes_;
}

const std::vector<NodeLabel>& Merge::GetResultNodeLabels() const {
  return result_node_labels_;
}

void Merge::ComputeResultEdgeMutations() {
  for (auto& [label, edge_id] : result_edges_) {
    Assert(label.GetParent().GetCompactGenome());
    Assert(label.GetChild().GetCompactGenome());
    const CompactGenome& parent = *label.GetParent().GetCompactGenome();
    const CompactGenome& child = *label.GetChild().GetCompactGenome();
    result_dag_.Get(edge_id).SetEdgeMutations(CompactGenome::ToEdgeMutations(
        result_dag_.GetReferenceSequence(), parent, child));
  }
}

bool Merge::ContainsLeafset(const LeafSet& leafset) const {
  return all_leaf_sets_.find(leafset) != all_leaf_sets_.end();
}

void Merge::ComputeLeafSets(const std::vector<size_t>& tree_idxs) {
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    MADAG tree = trees_.at(tree_idx);
    std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    std::vector<LeafSet> computed_ls = ComputeLeafSets(tree, labels);
    for (size_t node_idx = 0; node_idx < tree.GetNodesCount(); ++node_idx) {
      auto ls_iter = all_leaf_sets_.insert(std::move(computed_ls.at(node_idx)));
      labels.at(node_idx).SetLeafSet(std::addressof(*ls_iter.first));
    }
  });
}

void Merge::MergeTrees(const std::vector<size_t>& tree_idxs) {
  NodeId node_id{result_dag_.GetNodesCount()};
  std::mutex mtx;
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    for (auto label : labels) {
      std::unique_lock<std::mutex> lock{mtx};
      if (result_nodes_.try_emplace(label, node_id).second) {
        GetOrInsert(result_node_labels_, node_id) = label;
        ++node_id.value;
      }
    }
  });
  tbb::concurrent_vector<EdgeLabel> added_edges;
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    MADAG tree = trees_.at(tree_idx);
    const std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    for (Edge edge : tree.GetEdges()) {
      const auto& parent_label = labels.at(edge.GetParentId().value);
      const auto& child_label = labels.at(edge.GetChildId().value);
      auto ins = result_edges_.insert({{parent_label, child_label}, {}});
      if (ins.second) {
        added_edges.push_back({parent_label, child_label});
      }
    }
  });
  result_dag_.InitializeNodes(result_nodes_.size());
  EdgeId edge_id{result_dag_.GetEdgesCount()};
  for (auto edge : added_edges) {
    auto parent = result_nodes_.find(edge.GetParent());
    auto child = result_nodes_.find(edge.GetChild());
    Assert(parent != result_nodes_.end());
    Assert(child != result_nodes_.end());
    Assert(parent->second.value < result_dag_.GetNodesCount());
    Assert(child->second.value < result_dag_.GetNodesCount());
    result_dag_.AddEdge(edge_id, parent->second, child->second, edge.ComputeCladeIdx());
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

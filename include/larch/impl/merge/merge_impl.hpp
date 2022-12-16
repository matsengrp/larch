template <typename DAG>
Merge<DAG>::Merge(std::string_view reference_sequence) {
  ResultDAG().SetReferenceSequence(reference_sequence);
}

template <typename DAG>
void Merge<DAG>::AddDAGs(const std::vector<DAG>& dags) {
  for (auto dag : dags) {
    dag.AssertUA();
  }
  std::vector<size_t> dag_idxs;
  dag_idxs.resize(dags.size());
  std::iota(dag_idxs.begin(), dag_idxs.end(), trees_.size());

  for (DAG i : dags) {
    trees_.emplace_back(i);
  }
  tree_labels_.resize(trees_.size());

  tbb::parallel_for_each(dag_idxs.begin(), dag_idxs.end(), [&](size_t dag_idx) {
    DAG dag = trees_.at(dag_idx);
    std::vector<NodeLabel>& labels = tree_labels_.at(dag_idx);
    labels.resize(dag.GetNodesCount());
    for (size_t node_idx = 0; node_idx < dag.GetNodesCount(); ++node_idx) {
      auto cg_iter = ResultDAG().AddDeduplicated(
          dag.Get(NodeId{node_idx}).GetCompactGenome().Copy());
      labels.at(node_idx).SetCompactGenome(cg_iter.first);
    }
  });

  ComputeLeafSets(dag_idxs);
  MergeTrees(dag_idxs);
  Assert(result_nodes_.size() == ResultDAG().GetNodesCount());
  Assert(result_edges_.size() == ResultDAG().GetEdgesCount());
  ResultDAG().BuildConnections();
}

template <typename DAG>
template <typename D, typename N>
void Merge<DAG>::AddDAG(D dag, N below) {
  dag.AssertUA();
  tree_labels_.push_back({});
  std::vector<NodeLabel>& labels = tree_labels_.back();
  const bool is_subtree = not below.IsRoot();
  labels.resize(dag.GetNodesCount());
  for (auto node : dag.GetNodes()) {
    if (is_subtree and node.IsRoot()) {
      continue;
    }
    auto cg_iter = ResultDAG().AddDeduplicated(node.GetCompactGenome().Copy());
    labels.at(node.GetId().value).SetCompactGenome(cg_iter.first);
  }

  std::vector<LeafSet> computed_ls = ComputeLeafSets(dag, labels);
  for (auto node : dag.GetNodes()) {
    if (is_subtree and node.IsRoot()) {
      continue;
    }
    auto ls_iter = all_leaf_sets_.insert(std::move(computed_ls.at(node.GetId().value)));
    labels.at(node.GetId().value).SetLeafSet(std::addressof(*ls_iter.first));
  }

  NodeId node_id{ResultDAG().GetNodesCount()};
  for (auto label : labels) {
    if (is_subtree and label == NodeLabel{}) {
      continue;
    }
    if (result_nodes_.try_emplace(label, node_id).second) {
      GetOrInsert(result_node_labels_, node_id) = label;
      ++node_id.value;
    }
  }

  std::vector<EdgeLabel> added_edges;
  for (auto edge : dag.GetEdges()) {
    if (is_subtree and edge.IsRoot()) {
      continue;
    }
    const auto& parent_label = labels.at(edge.GetParentId().value);
    const auto& child_label = labels.at(edge.GetChildId().value);
    auto ins = result_edges_.insert({{parent_label, child_label}, {}});
    if (ins.second) {
      added_edges.push_back({parent_label, child_label});
    }
  }

  if (is_subtree) {
    EdgeLabel below_edge = {
        result_node_labels_.at(below.GetFirstParent().GetParent().GetId().value),
        labels.at(dag.GetRoot().GetFirstChild().GetChild().GetId().value)};
    auto ins = result_edges_.insert({below_edge, {}});
    if (ins.second) {
      added_edges.push_back(below_edge);
    }
  }

  ResultDAG().InitializeNodes(result_nodes_.size());
  EdgeId edge_id{ResultDAG().GetEdgesCount()};
  for (auto edge : added_edges) {
    auto parent = result_nodes_.find(edge.GetParent());
    auto child = result_nodes_.find(edge.GetChild());
    Assert(parent != result_nodes_.end());
    Assert(child != result_nodes_.end());
    Assert(parent->second.value < ResultDAG().GetNodesCount());
    Assert(child->second.value < ResultDAG().GetNodesCount());
    ResultDAG().AddEdge(edge_id, parent->second, child->second, edge.ComputeCladeIdx());
    ResultDAG().Get(parent->second) = parent->first.GetCompactGenome();
    ResultDAG().Get(child->second) = child->first.GetCompactGenome();
    auto result_edge_it = result_edges_.find(edge);
    Assert(result_edge_it != result_edges_.end());
    result_edge_it->second = edge_id;
    edge_id.value++;
  }

  Assert(result_nodes_.size() == ResultDAG().GetNodesCount());
  // TODO Assert(result_edges_.size() == ResultDAG().GetEdgesCount());
  ResultDAG().BuildConnections();
}

template <typename DAG>
MergeDAG Merge<DAG>::GetResult() const {
  return result_dag_storage_.View();
}

template <typename DAG>
MutableMergeDAG Merge<DAG>::ResultDAG() {
  return result_dag_storage_.View();
}

template <typename DAG>
const std::unordered_map<NodeLabel, NodeId>& Merge<DAG>::GetResultNodes() const {
  return result_nodes_;
}

template <typename DAG>
const std::vector<NodeLabel>& Merge<DAG>::GetResultNodeLabels() const {
  return result_node_labels_;
}

template <typename DAG>
void Merge<DAG>::ComputeResultEdgeMutations() {
  for (auto& [label, edge_id] : result_edges_) {
    Assert(label.GetParent().GetCompactGenome());
    Assert(label.GetChild().GetCompactGenome());
    const CompactGenome& parent = *label.GetParent().GetCompactGenome();
    const CompactGenome& child = *label.GetChild().GetCompactGenome();
    ResultDAG().Get(edge_id).SetEdgeMutations(CompactGenome::ToEdgeMutations(
        ResultDAG().GetReferenceSequence(), parent, child));
  }
}

template <typename DAG>
bool Merge<DAG>::ContainsLeafset(const LeafSet& leafset) const {
  return all_leaf_sets_.find(leafset) != all_leaf_sets_.end();
}

template <typename DAG>
void Merge<DAG>::ComputeLeafSets(const std::vector<size_t>& tree_idxs) {
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    DAG tree = trees_.at(tree_idx);
    std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    std::vector<LeafSet> computed_ls = ComputeLeafSets(tree, labels);
    for (size_t node_idx = 0; node_idx < tree.GetNodesCount(); ++node_idx) {
      auto ls_iter = all_leaf_sets_.insert(std::move(computed_ls.at(node_idx)));
      labels.at(node_idx).SetLeafSet(std::addressof(*ls_iter.first));
    }
  });
}

template <typename DAG>
void Merge<DAG>::MergeTrees(const std::vector<size_t>& tree_idxs) {
  NodeId node_id{ResultDAG().GetNodesCount()};
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
    DAG tree = trees_.at(tree_idx);
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
  ResultDAG().InitializeNodes(result_nodes_.size());
  EdgeId edge_id{ResultDAG().GetEdgesCount()};
  for (auto edge : added_edges) {
    auto parent = result_nodes_.find(edge.GetParent());
    auto child = result_nodes_.find(edge.GetChild());
    Assert(parent != result_nodes_.end());
    Assert(child != result_nodes_.end());
    Assert(parent->second.value < ResultDAG().GetNodesCount());
    Assert(child->second.value < ResultDAG().GetNodesCount());
    ResultDAG().AddEdge(edge_id, parent->second, child->second, edge.ComputeCladeIdx());
    ResultDAG().Get(parent->second) = parent->first.GetCompactGenome();
    ResultDAG().Get(child->second) = child->first.GetCompactGenome();
    auto result_edge_it = result_edges_.find(edge);
    Assert(result_edge_it != result_edges_.end());
    result_edge_it->second = edge_id;
    edge_id.value++;
  }
}

template <typename DAG>
std::vector<LeafSet> Merge<DAG>::ComputeLeafSets(DAG dag,
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

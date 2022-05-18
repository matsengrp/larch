Merge::Merge(const std::string& refseq,
             std::vector<std::reference_wrapper<const HistoryDAG>>&& trees,
             const std::vector<std::vector<NodeLabel>>& labels)
    : refseq_{refseq}, trees_{std::move(trees)}, labels_{labels} {
  assert(labels_.size() >= trees_.size());
  for (size_t tree_idx = 0; tree_idx < trees_.size(); ++tree_idx) {
    assert(labels_.at(tree_idx).size() >= trees_.at(tree_idx).get().GetNodes().size());
  }
}

void Merge::Run() {
  for (size_t tree_idx = 0; tree_idx < trees_.size(); ++tree_idx) {
    for (Edge edge : trees_.at(tree_idx).get().GetEdges()) {
      MakeResultEdge(tree_idx, edge.GetId());
    }
  }
  result_.BuildConnections();
}

const NodeLabel& Merge::GetNodeLabel(size_t tree_idx, NodeId node_id) {
  return labels_.at(tree_idx).at(node_id.value);
}

#include <iostream>

NodeId Merge::GetResultNode(const NodeLabel& label) {
  auto i = result_nodes_.find(label);
  NodeId id;
  if (i == result_nodes_.end()) {
    id = result_.AddNode({result_.GetNodes().size()}).GetId();
    result_nodes_.emplace_hint(i, label, id);
  } else {
    id = i->second;
  }
  return id;
}

void Merge::MakeResultEdge(size_t tree_idx, EdgeId edge_id) {
  Edge edge = trees_.at(tree_idx).get().GetEdge(edge_id);
  NodeId parent = edge.GetParent().GetId(), child = edge.GetChild().GetId();
  EdgeLabel label = {
      GetNodeLabel(tree_idx, parent), GetNodeLabel(tree_idx, child), {0}};
  auto i = result_edges_.find(label);
  if (i == result_edges_.end()) {
    NodeId result_parent = GetResultNode(std::get<0>(label));
    NodeId result_child = GetResultNode(std::get<1>(label));
    result_.AddEdge({result_.GetEdges().size()}, result_parent, result_child,
                    std::get<2>(label));
    result_edges_.emplace_hint(i, std::move(label));
  }
}

HistoryDAG& Merge::GetResult() { return result_; }

const HistoryDAG& Merge::GetResult() const { return result_; }

#include "mutation_annotated_dag.hpp"

const DAG& MADAG::GetDAG() const { return dag_; }

std::string_view MADAG::GetReferenceSequence() const { return reference_sequence_; }

const std::vector<EdgeMutations>& MADAG::GetEdgeMutations() const {
  return edge_mutations_;
}

const std::vector<CompactGenome>& MADAG::GetCompactGenomes() const {
  return compact_genomes_;
}

DAG& MADAG::GetDAG() { return dag_; }

std::string& MADAG::GetReferenceSequence() { return reference_sequence_; }

std::vector<EdgeMutations>& MADAG::GetEdgeMutations() { return edge_mutations_; }

std::vector<CompactGenome>& MADAG::GetCompactGenomes() { return compact_genomes_; }

std::vector<CompactGenome> MADAG::ComputeCompactGenomes(
    std::string_view reference_sequence_) const {
  Assert(not reference_sequence_.empty());
  Assert(not edge_mutations_.empty());
  std::vector<CompactGenome> result;
  result.resize(dag_.GetNodesCount());
  auto ComputeCG = [&](auto& self, Node node) {
    CompactGenome& compact_genome = result.at(node.GetId().value);
    if (node.IsRoot()) {
      compact_genome = CompactGenome(node, edge_mutations_, reference_sequence_);
      return;
    }
    if (not compact_genome.empty()) {
      return;
    }
    for (Edge edge : node.GetParents()) {
      self(self, edge.GetParent());
      const EdgeMutations& mutations = edge_mutations_.at(edge.GetId().value);
      const CompactGenome& parent = result.at(edge.GetParentId().value);
      compact_genome.AddParentEdge(mutations, parent, reference_sequence_);
    }
  };
  for (Node node : dag_.GetNodes()) {
    ComputeCG(ComputeCG, node);
  }
  return result;
}

std::vector<EdgeMutations> MADAG::ComputeEdgeMutations(
    std::string_view reference_sequence) const {
  Assert(not reference_sequence_.empty());
  Assert(not compact_genomes_.empty());
  std::vector<EdgeMutations> result;
  for (auto [parent, child] : dag_.GetEdges()) {
    result.emplace_back(CompactGenome::ToEdgeMutations(
        reference_sequence, compact_genomes_.at(parent.GetId().value),
        compact_genomes_.at(child.GetId().value)));
  }
  return result;
}

MADAG MADAG::GetSample() {
  MADAG sample_tree;
  sample_tree.GetReferenceSequence() = reference_sequence_;
  std::queue<Node> nodes_to_add;
  std::queue<Edge> edges_to_add;
  std::vector<NodeId> added_nodes;
  std::map<NodeId, size_t> new_node_idx;
  std::vector<EdgeId> old_edge_idx;
  size_t edge_ctr = 0;
  nodes_to_add.push(dag_.GetRoot());
  while (!nodes_to_add.empty()) {
    Node current_node = nodes_to_add.front();
    nodes_to_add.pop();
    added_nodes.emplace_back(current_node.GetId());
    for (auto clade : current_node.GetClades()) {
      int i = rand() % clade.size();
      Edge edge = clade[i];
      nodes_to_add.push(edge.GetChild());
      edges_to_add.push(edge);
    }
  }
  std::sort(added_nodes.begin(), added_nodes.end());
  for (size_t i = 0; i < added_nodes.size(); i++) {
    new_node_idx[added_nodes[i]] = i;
  }
  sample_tree.GetDAG().InitializeNodes(new_node_idx.size());
  while (!edges_to_add.empty()) {
    Edge edge = edges_to_add.front();
    edges_to_add.pop();
    sample_tree.GetDAG().AppendEdge({new_node_idx[edge.GetParent().GetId()]}, {new_node_idx[edge.GetChild().GetId()]}, edge.GetClade());
    old_edge_idx.emplace_back(edge.GetId());
  }
  sample_tree.GetDAG().BuildConnections();
  sample_tree.GetEdgeMutations().resize(old_edge_idx.size());
  for (size_t edge = 0; edge < old_edge_idx.size(); edge++) {
    for (auto mutation : MADAG::GetEdgeMutations().at(old_edge_idx[edge].value)) {
      sample_tree.GetEdgeMutations().at(edge)[mutation.first] = mutation.second;
    }
  }
  return sample_tree;
}

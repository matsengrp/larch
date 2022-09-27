#include "mutation_annotated_dag.hpp"
#include <string_view>

MADAG::MADAG(std::string_view reference_sequence)
    : reference_sequence_{reference_sequence} {}

const DAG& MADAG::GetDAG() const { return dag_; }

const std::string& MADAG::GetReferenceSequence() const { return reference_sequence_; }

const std::vector<EdgeMutations>& MADAG::GetEdgeMutations() const {
  return edge_mutations_;
}

const std::vector<CompactGenome>& MADAG::GetCompactGenomes() const {
  return compact_genomes_;
}

void MADAG::RemoveEdgeMutations() {
  edge_mutations_.resize(0);
  edge_mutations_.shrink_to_fit();
}

void MADAG::RemoveCompactGenomes() {
  compact_genomes_.resize(0);
  compact_genomes_.shrink_to_fit();
}

void MADAG::RecomputeCompactGenomes() { compact_genomes_ = ComputeCompactGenomes(); }

std::vector<CompactGenome> MADAG::ComputeCompactGenomes() const {
  return ComputeCompactGenomes(reference_sequence_, dag_, edge_mutations_);
}

std::vector<CompactGenome> MADAG::ComputeCompactGenomes(
    std::string_view reference_sequence, const DAG& dag,
    const std::vector<EdgeMutations>& edge_mutations) {
  Assert(not reference_sequence.empty());
  Assert(edge_mutations.size() == dag.GetEdgesCount());
  std::vector<CompactGenome> result;
  result.resize(dag.GetNodesCount());
  auto ComputeCG = [&](auto& self, Node node) {
    CompactGenome& compact_genome = result.at(node.GetId().value);
    if (node.IsRoot()) {
      compact_genome = CompactGenome(node, edge_mutations, reference_sequence);
      return;
    }
    if (not compact_genome.empty()) {
      return;
    }
    for (Edge edge : node.GetParents()) {
      self(self, edge.GetParent());
      const EdgeMutations& mutations = edge_mutations.at(edge.GetId().value);
      const CompactGenome& parent = result.at(edge.GetParentId().value);
      compact_genome.AddParentEdge(mutations, parent, reference_sequence);
    }
  };
  for (Node node : dag.GetNodes()) {
    ComputeCG(ComputeCG, node);
  }
  return result;
}

void MADAG::RecomputeEdgeMutations() { edge_mutations_ = ComputeEdgeMutations(); }

std::vector<EdgeMutations> MADAG::ComputeEdgeMutations() const {
  return ComputeEdgeMutations(reference_sequence_, dag_, compact_genomes_);
}

std::vector<EdgeMutations> MADAG::ComputeEdgeMutations(
    std::string_view reference_sequence, const DAG& dag,
    const std::vector<CompactGenome>& compact_genomes) {
  Assert(not reference_sequence.empty());
  Assert(compact_genomes.size() == dag.GetNodesCount());
  std::vector<EdgeMutations> result;
  for (auto [parent, child] : dag.GetEdges()) {
    result.emplace_back(CompactGenome::ToEdgeMutations(
        reference_sequence, compact_genomes.at(parent.GetId().value),
        compact_genomes.at(child.GetId().value)));
  }
  return result;
}

const EdgeMutations& MADAG::GetEdgeMutations(EdgeId edge_id) const {
  return edge_mutations_.at(edge_id.value);
}

MutableNode MADAG::AddNode(NodeId id) { return dag_.AddNode(id); }

MutableEdge MADAG::AddEdge(EdgeId id, NodeId parent, NodeId child, CladeIdx clade) {
  return dag_.AddEdge(id, parent, child, clade);
}

MutableEdge MADAG::AppendEdge(NodeId parent, NodeId child, CladeIdx clade) {
  return dag_.AppendEdge(parent, child, clade);
}

void MADAG::BuildConnections() { dag_.BuildConnections(); }

void MADAG::InitializeNodes(size_t nodes_count) { dag_.InitializeNodes(nodes_count); }

void MADAG::SetEdgeMutations(std::vector<EdgeMutations>&& edge_mutations) {
  MoveElements(std::forward<std::vector<EdgeMutations>>(edge_mutations),
               edge_mutations_);
}

void MADAG::AppendEdgeMutations(EdgeMutations&& edge_mutations) {
  edge_mutations_.push_back(std::forward<EdgeMutations>(edge_mutations));
}

void MADAG::AppendCompactGenome(CompactGenome&& compact_genome) {
  compact_genomes_.push_back(std::forward<CompactGenome>(compact_genome));
}

CompactGenome&& MADAG::ExtractCompactGenome(NodeId node) {
  return std::move(compact_genomes_.at(node.value));
}

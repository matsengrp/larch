#include "mutation_annotated_dag.hpp"

std::vector<CompactGenome> MADAG::ComputeCompactGenomes(
    std::string_view reference_sequence) const {
  std::vector<CompactGenome> result;
  result.resize(dag.GetNodes().size());
  for (auto [node, edge] : dag.TraversePreOrder()) {
    const EdgeMutations& mutations = edge_mutations.at(edge.GetId().value);
    const CompactGenome& parent = result.at(edge.GetParentId().value);
    CompactGenome& compact_genome = result.at(node.GetId().value);
    compact_genome = CompactGenome{mutations, parent, reference_sequence};
  }
  return result;
}

std::vector<CompactGenome> MADAG::ComputeCompactGenomesDAG(
    std::string_view reference_sequence) const {
  Assert(not reference_sequence.empty());
  Assert(not edge_mutations.empty());
  std::vector<CompactGenome> result;
  result.resize(dag.GetNodes().size());
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

std::vector<EdgeMutations> MADAG::ComputeEdgeMutations(
    std::string_view reference_sequence) const {
  Assert(not reference_sequence.empty());
  Assert(not compact_genomes.empty());
  std::vector<EdgeMutations> result;
  for (auto [parent, child] : dag.GetEdges()) {
    result.emplace_back(CompactGenome::ToEdgeMutations(
        reference_sequence, compact_genomes.at(parent.GetId().value),
        compact_genomes.at(child.GetId().value)));
  }
  return result;
}

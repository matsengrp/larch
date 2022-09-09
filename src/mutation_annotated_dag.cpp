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
      compact_genome = CompactGenome();
      return;
    }
    if (not compact_genome.empty()) {
      return;
    }
    Edge edge = *(node.GetParents().begin());
    self(self, edge.GetParent());
    const EdgeMutations& mutations = edge_mutations_.at(edge.GetId().value);
    const CompactGenome& parent = result.at(edge.GetParentId().value);
    compact_genome.AddParentEdge(mutations, parent, reference_sequence_);
  };
  std::unordered_map<CompactGenome, size_t> leaf_cgs;
  for (Node node : dag_.GetNodes()) {
    ComputeCG(ComputeCG, node);
    if (node.IsLeaf()) {
        bool success = leaf_cgs.emplace(result[node.GetId().value].Copy(), node.GetId().value).second;
        if (not success) {
            std::cout << "Error in ComputeCompactGenomes: had a non-unique leaf node at "
                << node.GetId().value
                << " also seen at "
                << leaf_cgs[result[node.GetId().value].Copy()]
                << "\nCompact Genome is\n"
                << result[node.GetId().value].ToString()
                << "\n" << std::flush;
            assert(false);
        }
    }

  }
  std::cout << "ComputeCompactGenomes found " << leaf_cgs.size() << " unique leaf cgs\n";
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

const EdgeMutations& MADAG::GetEdgeMutations(EdgeId edge_id) const {
  return edge_mutations_.at(edge_id.value);
}

EdgeMutations& MADAG::GetEdgeMutations(EdgeId edge_id) {
  return edge_mutations_.at(edge_id.value);
}

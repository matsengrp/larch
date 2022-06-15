#pragma once

#include <string_view>
#include <vector>

#include "dag.hpp"
#include "compact_genome.hpp"

class MADAG {
 public:
  const DAG& GetDAG() const;
  std::string_view GetReferenceSequence() const;
  const std::vector<EdgeMutations>& GetEdgeMutations() const;
  const std::vector<CompactGenome>& GetCompactGenomes() const;

  DAG& GetDAG();
  std::string& GetReferenceSequence();
  std::vector<EdgeMutations>& GetEdgeMutations();
  std::vector<CompactGenome>& GetCompactGenomes();

  [[nodiscard]] std::vector<CompactGenome> ComputeCompactGenomes(
      std::string_view reference_sequence) const;

  [[nodiscard]] std::vector<EdgeMutations> ComputeEdgeMutations(
      std::string_view reference_sequence) const;

 private:
  DAG dag_;
  std::string reference_sequence_;
  std::vector<EdgeMutations> edge_mutations_;
  std::vector<CompactGenome> compact_genomes_;
};

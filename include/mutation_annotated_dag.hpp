#pragma once

#include <string_view>
#include <vector>

#include "dag.hpp"
#include "compact_genome.hpp"

struct MADAG {
  DAG dag;
  std::string reference_sequence;
  std::vector<EdgeMutations> edge_mutations;
  std::vector<CompactGenome> compact_genomes;

  std::vector<CompactGenome> ComputeCompactGenomes(
      std::string_view reference_sequence) const;

  std::vector<CompactGenome> ComputeCompactGenomesDAG(
      std::string_view reference_sequence) const;

  std::vector<EdgeMutations> ComputeEdgeMutations(
      std::string_view reference_sequence) const;
};

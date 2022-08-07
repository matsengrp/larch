/**
 * MADAG extends DAG, containing additional node data (CompactGenomes) and edge
 * data (EdgeMutations) in vectors ordered by node/edge IDs.
 *
 * When edge_mutations_ is populated, MADAG implements the Mutation Annotated
 * DAG object, in which nucleotide sequences on nodes are determined by the
 * nucleotide sequence on the UA node (reference_sequence_) and a sequence of
 * mutations associated to any path of edges from the UA node to the node of
 * interest.
 *
 * In order to merge mutation annotated DAGs, nodes must be compared using
 * their sequence and child clade sets. Node sequences are represented by
 * CompactGenome objects, which may be computed from edge mutations using
 * GetCompactGenomes().
 *
 */

#pragma once

#include <string_view>
#include <vector>
#include <queue>

#include "dag.hpp"
#include "compact_genome.hpp"

class MADAG {
 public:
  const DAG& GetDAG() const;
  /**
   * Get UA node sequence
   */
  std::string_view GetReferenceSequence() const;
  /**
   * Edge mutations are ordered by edge ID
   */
  const std::vector<EdgeMutations>& GetEdgeMutations() const;
  /**
   * Compact genomes are ordered by node ID
   */
  const std::vector<CompactGenome>& GetCompactGenomes() const;

  DAG& GetDAG();
  std::string& GetReferenceSequence();
  std::vector<EdgeMutations>& GetEdgeMutations();
  std::vector<CompactGenome>& GetCompactGenomes();

  MADAG GetSample();
  /**
   * Compute compact genomes, ordered by node ID, from edge mutations
   */
  [[nodiscard]] std::vector<CompactGenome> ComputeCompactGenomes(
      std::string_view reference_sequence) const;

  /**
   * Compute edge mutations, ordered by edge ID, from compact genomes
   */
  [[nodiscard]] std::vector<EdgeMutations> ComputeEdgeMutations(
      std::string_view reference_sequence) const;

 private:
  DAG dag_;
  std::string reference_sequence_;
  std::vector<EdgeMutations> edge_mutations_;
  std::vector<CompactGenome> compact_genomes_;
};

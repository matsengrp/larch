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
#include <unordered_map>
#include <utility>

#include "dag.hpp"
#include "compact_genome.hpp"
#include "edge.hpp"

namespace Mutation_Annotated_Tree {
class Tree;
class Node;
}  // namespace Mutation_Annotated_Tree

class MADAG {
 public:
  MADAG() = default;
  explicit MADAG(std::string_view reference_sequence);

  const DAG& GetDAG() const;
  DAG& GetDAG();
  /**
   * Get UA node sequence
   */
  const std::string& GetReferenceSequence() const;
  /**
   * Edge mutations are ordered by edge ID
   */
  const std::vector<EdgeMutations>& GetEdgeMutations() const;
  /**
   * Compact genomes are ordered by node ID
   */
  const std::vector<CompactGenome>& GetCompactGenomes() const;

  /**
   * Clears all stored edge mutations
   */
  void RemoveEdgeMutations();

  /**
   * Clears all stored compact genomes
   */
  void RemoveCompactGenomes();

  /**
   * Compute compact genomes from stored edge mutatons, and store internally
   */
  void RecomputeCompactGenomes();

  /**
   * Compute compact genomes from stored edge mutatons, doesn't change the existing
   * stored mutations
   */
  [[nodiscard]] std::vector<CompactGenome> ComputeCompactGenomes() const;

  /**
   * Compute compact genomes, ordered by node ID, from edge mutations
   */
  [[nodiscard]] static std::vector<CompactGenome> ComputeCompactGenomes(
      std::string_view reference_sequence, const DAG& dag,
      const std::vector<EdgeMutations>& edge_mutations);
  /**
   * Compute edge mutations from stored compact genomes, and store internally
   */
  void RecomputeEdgeMutations();

  /**
   * Compute edge mutations from stored compact genomes, doesn't change the existing
   * stored mutations
   */
  [[nodiscard]] std::vector<EdgeMutations> ComputeEdgeMutations() const;

  /**
   * Compute edge mutations, ordered by edge ID, from compact genomes
   */
  [[nodiscard]] static std::vector<EdgeMutations> ComputeEdgeMutations(
      std::string_view reference_sequence, const DAG& dag,
      const std::vector<CompactGenome>& compact_genomes);

  const EdgeMutations& GetEdgeMutations(EdgeId edge_id) const;

  bool HaveUA() const;
  void AssertUA() const;
  void AddUA(const EdgeMutations& mutations_at_root);

 private:
  friend MADAG LoadDAGFromProtobuf(std::string_view);
  friend MADAG LoadTreeFromProtobuf(std::string_view, std::string_view);
  friend MADAG LoadDAGFromJson(std::string_view);
  friend MADAG MakeSyntheticDAG();
  template <typename>
  friend class SubtreeWeight;
  friend class Merge;
  friend MADAG build_madag_from_mat(const Mutation_Annotated_Tree::Tree&,
                                    std::string_view);
  friend void build_madag_from_mat_helper(Mutation_Annotated_Tree::Node*, Node, MADAG&);

  MutableNode AddNode(NodeId id);
  MutableNode AppendNode();
  MutableEdge AddEdge(EdgeId id, NodeId parent, NodeId child, CladeIdx clade);
  MutableEdge AppendEdge(NodeId parent, NodeId child, CladeIdx clade);
  void BuildConnections();
  void InitializeNodes(size_t nodes_count);
  void SetEdgeMutations(std::vector<EdgeMutations>&& edge_mutations);
  void AppendEdgeMutations(EdgeMutations&& edge_mutations);
  void AppendCompactGenome(CompactGenome&& compact_genome);
  CompactGenome&& ExtractCompactGenome(NodeId node);
  void ResizeEdgeMutations(size_t size);
  void SetEdgeMutations(EdgeId edge, EdgeMutations&& mutations);

  DAG dag_;
  std::string reference_sequence_;
  std::vector<EdgeMutations> edge_mutations_;
  std::vector<CompactGenome> compact_genomes_;
};

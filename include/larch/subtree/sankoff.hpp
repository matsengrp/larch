/**
 * Sankoff algorithm implementation for weighted parsimony scoring.
 *
 * Supports arbitrary per-site 4x4 substitution cost matrices.
 * Computes weighted parsimony scores via bottom-up DP and
 * reconstructs ancestral sequences via top-down traceback.
 */

#pragma once

#include <array>
#include <limits>
#include <unordered_map>
#include <vector>

#include "larch/madag/mutation_annotated_dag.hpp"

// 4x4 substitution cost matrix for a single site
// Indexed by [parent_base][child_base] where base indices are A=0, C=1, G=2, T=3
using SiteCostMatrix = std::array<std::array<double, 4>, 4>;

// Cost array for a single site at a node: cost[base] = min cost to explain subtree
// if this node has the given base
using SiteCosts = std::array<double, 4>;

// Infinity value for impossible states
constexpr double kSankoffInfinity = std::numeric_limits<double>::infinity();

// Convert base character to index (A=0, C=1, G=2, T=3)
inline size_t BaseToIndex(char base);

// Convert index to base character
inline char IndexToBase(size_t index);

// Get indices of all bases compatible with a MutationBase (handles ambiguity)
inline std::vector<size_t> GetCompatibleIndices(MutationBase base);

/**
 * DP table for Sankoff algorithm.
 * Stores costs and traceback information for variable sites only.
 */
struct SankoffDPTable {
  // Variable sites that need scoring (positions that have mutations)
  std::vector<MutationPosition> variable_sites;

  // Map from position to index in variable_sites
  std::unordered_map<size_t, size_t> site_to_index;

  // DP costs: dp_costs[node_id][site_index] = SiteCosts (cost for each base)
  std::vector<std::vector<SiteCosts>> dp_costs;

  // Traceback: traceback[edge_id][site_index][parent_base] = optimal child base index
  std::vector<std::vector<std::array<uint8_t, 4>>> traceback;

  // Per-site cost matrices (indexed by site_index)
  std::vector<SiteCostMatrix> cost_matrices;
};

/**
 * Sankoff algorithm scorer for weighted parsimony on trees.
 *
 * @tparam DAG The DAG type to score (must support MADAG interface)
 *
 * Usage:
 *   auto scorer = SankoffScorer(dag);
 *   double score = scorer.ComputeScoreBelow(dag.GetRoot());
 *   scorer.ReconstructAncestralSequences(dag.GetRoot());
 *   auto genome = scorer.GetReconstructedGenome(node);
 */
template <typename DAG>
class SankoffScorer {
 public:
  /**
   * Construct scorer with uniform costs (0 diagonal, 1 off-diagonal).
   */
  explicit SankoffScorer(DAG dag);

  /**
   * Construct scorer with site-specific cost matrices.
   * @param cost_matrices Map from MutationPosition to cost matrix. Sites not in
   *                      the map use uniform costs.
   */
  SankoffScorer(DAG dag,
                const std::unordered_map<size_t, SiteCostMatrix>& cost_matrices);

  /**
   * Construct scorer with a single cost matrix applied to all sites.
   */
  SankoffScorer(DAG dag, const SiteCostMatrix& uniform_matrix);

  /**
   * Collect all variable sites from edge mutations.
   * Called automatically by constructor.
   */
  void CollectVariableSites();

  /**
   * Compute weighted parsimony score below the given node (bottom-up DP).
   * @return Total score (sum of minimum costs across all variable sites at root)
   */
  double ComputeScoreBelow(typename DAG::NodeView node);

  /**
   * Reconstruct ancestral sequences via top-down traceback.
   * Must be called after ComputeScoreBelow().
   */
  void ReconstructAncestralSequences(typename DAG::NodeView root);

  /**
   * Get the reconstructed genome for a node as CompactGenome.
   * Returns mutations relative to the reference sequence.
   * Must be called after ReconstructAncestralSequences().
   */
  CompactGenome GetReconstructedGenome(typename DAG::NodeView node) const;

  /**
   * Get the reconstructed base at a specific site for a node.
   * @param site_index Index into variable_sites (not the position value)
   */
  char GetReconstructedBase(NodeId node_id, size_t site_index) const;

  /**
   * Get total score (cached from last ComputeScoreBelow call).
   */
  double GetTotalScore() const { return total_score_; }

  /**
   * Get number of variable sites being tracked.
   */
  size_t GetNumVariableSites() const { return dp_table_.variable_sites.size(); }

  /**
   * Get const reference to DP table (for debugging/inspection).
   */
  const SankoffDPTable& GetDPTable() const { return dp_table_; }

 private:
  DAG dag_;
  SankoffDPTable dp_table_;
  double total_score_ = 0.0;
  bool score_computed_ = false;
  bool ancestors_reconstructed_ = false;

  // Reconstructed bases: ancestral_bases_[node_id][site_index] = base index (0-3)
  std::vector<std::vector<uint8_t>> ancestral_bases_;

  // Per-site cost matrices provided at construction (position -> matrix)
  std::unordered_map<size_t, SiteCostMatrix> custom_cost_matrices_;

  // Single uniform matrix to use for all sites (if provided)
  std::optional<SiteCostMatrix> uniform_cost_matrix_;

  /**
   * Initialize cost matrices for all variable sites.
   */
  void InitializeCostMatrices();

  /**
   * Create default uniform cost matrix (0 diagonal, 1 off-diagonal).
   */
  static SiteCostMatrix MakeUniformCostMatrix();

  /**
   * Compute DP costs for a single node (recursive, bottom-up).
   */
  void ComputeNodeCosts(typename DAG::NodeView node);

  /**
   * Traceback from a node to assign bases to descendants.
   * @param assigned_bases Base indices assigned to this node for each site
   */
  void TracebackFromNode(typename DAG::NodeView node,
                         const std::vector<uint8_t>& assigned_bases);

  /**
   * Get the observed base at a leaf node for a given site.
   * Returns the base from CompactGenome or reference sequence.
   */
  MutationBase GetLeafBase(typename DAG::NodeView leaf, MutationPosition pos) const;
};

#include "larch/impl/subtree/sankoff_impl.hpp"

#pragma once

#include "larch/spr/score.hpp"
#include "larch/contiguous_set.hpp"
#include "larch/contiguous_map.hpp"
#include "larch/usher_glue.hpp"

#ifdef USE_NETAM
#include <netam/crepe.hpp>
#include <netam/likelihood.hpp>
#include <netam/kmer_sequence_encoder.hpp>
#endif

/**
 * @brief FitchSet represents a set of possible nucleotides at a site (for Fitch
 * algorithm).
 *
 * Internally uses a 4-bit one-hot encoding where:
 * - bit 0 (1) = A
 * - bit 1 (2) = C
 * - bit 2 (4) = G
 * - bit 3 (8) = T
 */
struct FitchSet {
  FitchSet(char base) : value_{base} { Assert(value_ > 0); }
  FitchSet(int base) : value_{static_cast<char>(base)} { Assert(value_ > 0); }
  bool find(char base) const {
    switch (base) {
      case 'A':
        return (value_ & 1) != 0;
      case 'C':
        return (value_ & 2) != 0;
      case 'G':
        return (value_ & 4) != 0;
      case 'T':
        return (value_ & 8) != 0;  // NOLINT
    }
    Fail("Unreachable");
  }
  char at([[maybe_unused]] size_t pos) const {
    Assert(pos == 0);
    if ((value_ & 1) != 0) {
      return 'A';
    } else if ((value_ & 2) != 0) {
      return 'C';
    } else if ((value_ & 4) != 0) {
      return 'G';
    } else if ((value_ & 8) != 0) {  // NOLINT
      return 'T';
    }
    Fail("Unreachable");
  }

 private:
  char value_ = 0;
};

inline nuc_one_hot base_to_singleton(MutationBase base) {
  switch (base.ToChar()) {
    case 'A':
      return 1;
    case 'C':
      return 2;
    case 'G':
      return 4;
    case 'T':
      return 8;  // NOLINT
  }
  Fail("unrecognized base");
}

/**
 * @brief MatOptimize-based scoring backend using Fitch algorithm.
 *
 * This backend uses the existing matOptimize infrastructure to compute
 * parsimony scores based on Fitch sets and mutation count changes.
 *
 * The backend stores:
 * - Move information (src, dst, lca as NodeId)
 * - Changed Fitch sets keyed by NodeId (converted from MATNodePtr)
 * - Reference to the DAG for MAT node lookups
 */
template <typename DAG>
class MatOptimizeScoringBackend {
 public:
  using NucleotideSet = FitchSet;

  MatOptimizeScoringBackend() = default;

  /**
   * @brief Initialize from matOptimize move data.
   *
   * Converts MAT::Node* to NodeId and stores scoring data internally.
   * Must be called after ApplyMove() has been applied to the DAG.
   */
  template <typename DAGView>
  bool Initialize(const DAGView& dag, const Profitable_Moves& move,
                  const std::vector<Node_With_Major_Allele_Set_Change>& changes);

  // Move access (always returns NodeId)
  NodeId GetMoveSource() const { return src_; }
  NodeId GetMoveTarget() const { return dst_; }
  NodeId GetMoveLCA() const { return lca_; }
  Score GetScoreChange() const { return Score{move_score_change_}; }

  // Per-node scoring queries
  ContiguousSet<MutationPosition> GetSitesWithScoringChanges(NodeId node) const;
  bool HasScoringChanges(NodeId node) const;

  /**
   * @brief Get the scoring data for a node.
   *
   * Returns the Fitch set parts (mutations collection and optional changes map)
   * for the given node, accounting for move-related node transformations.
   */
  template <typename DAGView>
  std::pair<MAT::Mutations_Collection,
            std::optional<ContiguousMap<MutationPosition, Mutation_Count_Change>>>
  GetFitchSetParts(const DAGView& dag, NodeId node, bool is_leaf, bool is_move_target,
                   bool is_move_new) const;

  /**
   * @brief Compute the Fitch set at a specific site for a node.
   */
  template <typename DAGView>
  FitchSet GetFitchSetAtSite(const DAGView& dag, NodeId node, MutationPosition site,
                             bool is_leaf, bool is_move_target, bool is_move_new) const;

  /**
   * @brief Select the best base from a Fitch set given constraints.
   */
  static MutationBase SelectBase(const FitchSet& fitch_set, MutationBase old_base,
                                 MutationBase parent_base);

  // Access to internal data for legacy compatibility
  const ContiguousMap<NodeId, ContiguousMap<MutationPosition, Mutation_Count_Change>>&
  GetChangedFitchSetMap() const {
    return changed_fitch_set_map_;
  }

 private:
  // Converted move info (NodeId-based)
  NodeId src_;
  NodeId dst_;
  NodeId lca_;
  int move_score_change_ = 0;

  // Changed Fitch sets, keyed by NodeId (converted from MATNodePtr during init)
  ContiguousMap<NodeId, ContiguousMap<MutationPosition, Mutation_Count_Change>>
      changed_fitch_set_map_;

  // Keep original MATNodePtr-keyed map for access through MAT
  ContiguousMap<MATNodePtr, ContiguousMap<MutationPosition, Mutation_Count_Change>>
      mat_keyed_fitch_set_map_;
};

/**
 * @brief ML-based scoring backend using neural network likelihood.
 *
 * This backend uses netam::crepe (which wraps indep_rscnn_model and
 * kmer_sequence_encoder) to compute log-likelihoods for edges in the
 * affected subtree of an SPR move.
 *
 * Key design decisions:
 * - Edge-level likelihood only (sum log-likelihoods for affected edges)
 * - Keep using parsimony (Fitch) for internal node base selection
 * - Score is negated log-likelihood (lower is better, consistent with parsimony)
 * - Simple implementation without batching
 */
template <typename DAG>
class MLScoringBackend {
 public:
  /**
   * @brief Simple nucleotide set for ML backend - just holds a reference base.
   * Base selection still uses parsimony logic.
   */
  struct NucleotideSet {
    char base;
    bool find(char b) const { return base == b; }
    char at(size_t) const { return base; }
  };

  MLScoringBackend() = default;

#ifdef USE_NETAM
  /**
   * @brief Construct with netam model for ML scoring.
   */
  explicit MLScoringBackend(netam::crepe& model) : model_(&model) {}
#endif

  /**
   * @brief Initialize from NodeId-based move and compute ML score.
   */
  template <typename DAGView>
  bool Initialize(const DAGView& dag, NodeId src, NodeId dst, NodeId lca);

  // Move access
  NodeId GetMoveSource() const { return src_; }
  NodeId GetMoveTarget() const { return dst_; }
  NodeId GetMoveLCA() const { return lca_; }

  /**
   * @brief Get parsimony score change (always 0 for ML backend).
   * Use GetMLScoreChange() for the actual ML score.
   */
  Score GetScoreChange() const { return Score{0}; }

  /**
   * @brief Get ML score change (negated log-likelihood difference).
   * Lower is better, consistent with parsimony scoring.
   */
  double GetMLScoreChange() const { return ml_score_change_; }

  // Per-node scoring queries - empty for ML backend (uses parsimony for base selection)
  ContiguousSet<MutationPosition> GetSitesWithScoringChanges(NodeId) const {
    return {};
  }
  bool HasScoringChanges(NodeId) const { return false; }

  /**
   * @brief Get scoring data for a node (returns empty - ML doesn't use Fitch sets).
   */
  template <typename DAGView>
  std::pair<MAT::Mutations_Collection,
            std::optional<ContiguousMap<MutationPosition, Mutation_Count_Change>>>
  GetFitchSetParts(const DAGView&, NodeId, bool, bool, bool) const {
    return {MAT::Mutations_Collection{}, std::nullopt};
  }

  /**
   * @brief Get nucleotide set at site (returns base from compact genome).
   */
  template <typename DAGView>
  NucleotideSet GetNucleotideSetAtSite(const DAGView& dag, NodeId node,
                                       MutationPosition site, bool, bool, bool) const;

  /**
   * @brief Select base (keeps old base - parsimony handles base selection).
   */
  static MutationBase SelectBase(const NucleotideSet&, MutationBase old_base,
                                 MutationBase) {
    return old_base;
  }

#ifdef USE_NETAM
  /**
   * @brief Expand CompactGenome to full sequence string.
   */
  template <typename DAGView>
  static std::string ExpandSequence(const DAGView& dag, NodeId node);

  /**
   * @brief Compute log-likelihood for a single edge (parent -> child).
   */
  template <typename DAGView>
  double ComputeEdgeLogLikelihood(const DAGView& dag, NodeId parent,
                                  NodeId child) const;
#endif

 private:
  NodeId src_;
  NodeId dst_;
  NodeId lca_;
  double ml_score_change_ = 0.0;

#ifdef USE_NETAM
  netam::crepe* model_ = nullptr;
#endif
};

/**
 * @brief Pure-MADAG scoring backend using Fitch algorithm.
 *
 * This backend computes Fitch sets via bottom-up parsimony directly on the
 * overlay tree created by ApplyMove. No MAT dependency — works entirely
 * with larch-native NodeIds and DAG traversal.
 *
 * The backend runs a full Fitch (unweighted parsimony) pass on both the
 * pre-move and post-move trees during Initialize, comparing the results
 * to identify which sites changed their Fitch set at each node.
 */
template <typename DAG>
class ParsimonyOnlyScoringBackend {
 public:
  using NucleotideSet = FitchSet;

  ParsimonyOnlyScoringBackend() = default;

  /**
   * @brief Initialize by running Fitch algorithm on overlay and original trees.
   *
   * Computes Fitch sets for all variable sites at all nodes, compares
   * pre-move and post-move results to identify changes and score delta.
   */
  template <typename DAGView>
  bool Initialize(const DAGView& dag, NodeId src, NodeId dst, NodeId lca);

  // Move access (always returns NodeId)
  NodeId GetMoveSource() const { return src_; }
  NodeId GetMoveTarget() const { return dst_; }
  NodeId GetMoveLCA() const { return lca_; }
  Score GetScoreChange() const { return Score{score_change_}; }

  // Per-node scoring queries
  ContiguousSet<MutationPosition> GetSitesWithScoringChanges(NodeId node) const;
  bool HasScoringChanges(NodeId node) const;

  /**
   * @brief Get scoring data for a node.
   *
   * Returns empty Mutations_Collection and a changes map whose keys are the
   * sites where the Fitch set changed. Values are default Mutation_Count_Change
   * (only the keys matter for ComputeNewCompactGenome).
   */
  template <typename DAGView>
  std::pair<MAT::Mutations_Collection,
            std::optional<ContiguousMap<MutationPosition, Mutation_Count_Change>>>
  GetFitchSetParts(const DAGView& dag, NodeId node, bool is_leaf, bool is_move_target,
                   bool is_move_new) const;

  /**
   * @brief Compute the Fitch set at a specific site for a node.
   *
   * Returns the pre-computed Fitch set from the bottom-up pass.
   */
  template <typename DAGView>
  FitchSet GetFitchSetAtSite(const DAGView& dag, NodeId node, MutationPosition site,
                             bool is_leaf, bool is_move_target, bool is_move_new) const;

  /**
   * @brief Select the best base from a Fitch set given constraints.
   * Same logic as MatOptimizeScoringBackend.
   */
  static MutationBase SelectBase(const FitchSet& fitch_set, MutationBase old_base,
                                 MutationBase parent_base);

  // Returns empty map (no MAT data)
  const ContiguousMap<MATNodePtr, ContiguousMap<MutationPosition, Mutation_Count_Change>>&
  GetChangedFitchSetMap() const;

 private:
  NodeId src_;
  NodeId dst_;
  NodeId lca_;
  int score_change_ = 0;

  // Ordered variable sites (positions with mutations in the tree)
  std::vector<MutationPosition> variable_sites_;

  // Map from position value to index in variable_sites_
  std::unordered_map<size_t, size_t> site_to_index_;

  // Post-move Fitch sets: [node_id.value][site_index] = nuc_one_hot bitmask
  std::unordered_map<size_t, std::vector<uint8_t>> new_fitch_sets_;

  // Sites with changed Fitch sets per node
  ContiguousMap<NodeId, ContiguousSet<MutationPosition>> changed_sites_map_;

  // Recursive Fitch bottom-up pass
  template <typename TreeView>
  void FitchVisit(const TreeView& tree, NodeId node_id,
                  std::unordered_map<size_t, std::vector<uint8_t>>& fitch_sets,
                  int& score) const;

  // Find the tree root (non-UA node whose parent is UA)
  template <typename TreeView>
  NodeId FindTreeRoot(const TreeView& tree) const;
};

// Implementation included at the end
#include "larch/impl/spr/scoring_backend_impl.hpp"

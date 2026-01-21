#pragma once

#include "larch/spr/score.hpp"
#include "larch/contiguous_set.hpp"
#include "larch/contiguous_map.hpp"
#include "larch/usher_glue.hpp"

/**
 * @brief FitchSet represents a set of possible nucleotides at a site (for Fitch algorithm).
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
 * @brief ML/Sankoff-based scoring backend (stub implementation).
 *
 * This backend is a placeholder for future neural network-based likelihood
 * scoring. Currently returns Score{0} for all operations.
 *
 * The backend stores only move topology (src, dst, lca as NodeId) since
 * scoring is deferred to the ML model.
 */
template <typename DAG>
class MLScoringBackend {
 public:
  /**
   * @brief Simple nucleotide set for ML backend - just holds a reference base.
   */
  struct NucleotideSet {
    char base;
    bool find(char b) const { return base == b; }
    char at(size_t) const { return base; }
  };

  MLScoringBackend() = default;

  /**
   * @brief Initialize from NodeId-based move.
   */
  template <typename DAGView>
  bool Initialize(const DAGView& dag, NodeId src, NodeId dst, NodeId lca);

  // Move access
  NodeId GetMoveSource() const { return src_; }
  NodeId GetMoveTarget() const { return dst_; }
  NodeId GetMoveLCA() const { return lca_; }
  Score GetScoreChange() const { return Score{0}; }  // Stub: always 0

  // Per-node scoring queries - empty for stub
  ContiguousSet<MutationPosition> GetSitesWithScoringChanges(NodeId) const { return {}; }
  bool HasScoringChanges(NodeId) const { return false; }

  /**
   * @brief Get scoring data for a node (stub - returns empty/default).
   */
  template <typename DAGView>
  std::pair<MAT::Mutations_Collection,
            std::optional<ContiguousMap<MutationPosition, Mutation_Count_Change>>>
  GetFitchSetParts(const DAGView&, NodeId, bool, bool, bool) const {
    return {MAT::Mutations_Collection{}, std::nullopt};
  }

  /**
   * @brief Get nucleotide set at site (stub - returns reference base).
   */
  template <typename DAGView>
  NucleotideSet GetNucleotideSetAtSite(const DAGView& dag, NodeId node, MutationPosition site,
                                       bool, bool, bool) const;

  /**
   * @brief Select base (stub - keeps old base).
   */
  static MutationBase SelectBase(const NucleotideSet&, MutationBase old_base,
                                 MutationBase) {
    return old_base;
  }

 private:
  NodeId src_;
  NodeId dst_;
  NodeId lca_;
};

// Implementation included at the end
#include "larch/impl/spr/scoring_backend_impl.hpp"

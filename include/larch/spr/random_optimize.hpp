#pragma once

#include <iostream>
#include <optional>
#include <vector>

#include "larch/merge/merge.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/spr/spr_view.hpp"
#include "larch/spr/scoring_backend.hpp"
#include "larch/spr/random_moves.hpp"

/**
 * @brief Pure-MADAG optimization loop using random SPR moves.
 *
 * Implements the sample->move->fragment->merge loop using only larch-native
 * MADAG structures and NodeIds. No matOptimize dependency.
 *
 * Each iteration:
 * 1. Sample a tree from the DAG
 * 2. Generate random SPR moves on the sampled tree
 * 3. For each move, create a hypothetical tree and extract a fragment
 * 4. Merge all fragments (and the sampled tree) into the DAG
 *
 * @param merge The Merge object to grow
 * @param num_iterations Number of sample-optimize-merge iterations
 * @param moves_per_iteration Number of random moves to attempt per iteration
 * @param collapse_empty_fragment_edges Whether to collapse edges with no mutations
 * @param seed Optional random seed for reproducibility
 */
inline void OptimizeDAGWithRandomMoves(
    Merge& merge, size_t num_iterations, size_t moves_per_iteration,
    bool collapse_empty_fragment_edges = true,
    std::optional<uint32_t> seed = std::nullopt) {
  std::ignore = collapse_empty_fragment_edges;

  std::mt19937 seed_gen;
  if (seed.has_value()) {
    seed_gen.seed(seed.value());
  } else {
    std::random_device rd;
    seed_gen.seed(rd());
  }

  for (size_t iter = 0; iter < num_iterations; iter++) {
    // 1. Compute edge mutations on the merged result
    merge.ComputeResultEdgeMutations();

    // 2. Sample a min-weight (parsimony) tree from the DAG
    SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
    auto sampled_storage = weight.MinWeightSampleTree({});
    auto sampled = sampled_storage.View();
    sampled.RecomputeCompactGenomes(true);
    sampled.SampleIdsFromCG();

    // Verify sampled tree is valid
    Assert(sampled.IsTree());

    // 3. Generate random moves on the sampled tree
    uint32_t move_seed = seed_gen();
    using SampledView = decltype(sampled.Const());
    RandomMoveGenerator<SampledView> move_gen{sampled.Const(), move_seed};

    using SampledStorage = std::remove_reference_t<decltype(sampled_storage)>;
    using Backend = ParsimonyOnlyScoringBackend<SampledStorage>;

    size_t accepted_moves = 0;

    for (size_t m = 0; m < moves_per_iteration; m++) {
      auto move = move_gen.GenerateMove();
      if (not move.has_value()) {
        continue;
      }

      // Create SPR storage with ParsimonyOnly backend (mutable view needed
      // for overlay features like Deduplicate<CompactGenome>)
      auto spr = AddSPRStorageWithBackend<Backend>(sampled);

      // Apply move and initialize hypothetical tree
      bool success =
          spr.View().InitHypotheticalTree(move->src, move->dst, move->lca);
      if (not success) {
        continue;
      }

      // Extract fragment and overlay MappedNodes on SPR nodes so that
      // merge can write SetOriginalId through the fragment view.
      // (SampledDAGStorage has MappedNodes, which propagates through the SPR
      // overlay to the fragment. Without overlaying, writes to MappedNodes on
      // non-overlaid nodes fail with "Can't modify non-overlaid node".)
      auto fragment = spr.View().MakeFragment();
      for (auto node : fragment.View().GetNodes()) {
        auto spr_node = spr.View().Get(node.GetId());
        if (not spr_node.IsAppended()) {
          spr_node.template SetOverlay<MappedNodes>();
        }
      }
      merge.AddDAG(fragment.View());
      accepted_moves++;
    }

    // 4. Merge the sampled tree itself (mutable view needed for MappedNodes)
    merge.AddDAG(sampled);
  }
}

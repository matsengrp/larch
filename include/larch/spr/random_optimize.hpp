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
#include "larch/parallel/reduction.hpp"
#include "larch/benchmark.hpp"

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

    [[maybe_unused]] size_t accepted_moves = 0;

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

/**
 * @brief Compute max depth of a tree (max distance from tree root to any leaf).
 */
template <typename DAG>
size_t ComputeTreeMaxDepth(DAG dag) {
  size_t max_depth = 0;
  for (auto node : dag.GetNodes()) {
    if (not node.IsLeaf()) {
      continue;
    }
    size_t depth = 0;
    auto current = node;
    while (not current.IsUA()) {
      depth++;
      current = current.GetSingleParent().GetParent();
    }
    if (depth > max_depth) {
      max_depth = depth;
    }
  }
  return max_depth;
}

/**
 * @brief Result of a single parallel optimization radius pass.
 */
struct RadiusResult {
  size_t accepted_moves;
  long apply_ms;     // Phase 1 (PARALLEL): move generation + SPR overlay
  long fragment_ms;  // Phase 2 (PARALLEL): MakeFragment + SetOverlay<MappedNodes>
  long merge_ms;     // Phase 3 (BATCH):    AddDAGs (internally parallel)
};

/**
 * @brief Parallel MADAG-only optimization using random SPR moves with radius iteration.
 *
 * Similar to the matOptimize-based optimize_dag_direct, but operates entirely on
 * MADAG structures. Three-phase pipeline per radius:
 *   Phase 1 (PARALLEL): Generate moves + apply SPR overlays
 *   Phase 2 (PARALLEL): Create fragments + overlay MappedNodes
 *   Phase 3 (BATCH):    Merge all fragments via AddDAGs (internally parallel)
 *
 * @param merge The Merge object to grow
 * @param sampled_storage The sampled tree storage (will be used for move generation)
 * @param moves_per_radius Number of random moves to attempt per radius
 * @param seed Random seed for reproducibility
 * @return Vector of RadiusResult for timing analysis
 */
template <typename SampledStorage>
std::vector<RadiusResult> OptimizeDAGParallelRadius(
    Merge& merge, SampledStorage& sampled_storage, size_t moves_per_radius,
    uint32_t seed) {
  auto sampled = sampled_storage.View();
  auto sampled_const = sampled.Const();

  using SampledConstView = decltype(sampled_const);
  using Backend = ParsimonyOnlyScoringBackend<std::remove_reference_t<SampledStorage>>;

  // Compute max radius from tree depth
  size_t max_depth = ComputeTreeMaxDepth(sampled_const);
  size_t max_radius = max_depth * 2;
  std::cout << "maximum radius is " << max_radius << "\n" << std::flush;

  std::vector<RadiusResult> results;

  for (size_t rad_exp = 1; (static_cast<size_t>(1) << rad_exp) <= max_radius;
       rad_exp++) {
    size_t radius = static_cast<size_t>(1) << rad_exp;
    std::cout << "current radius is " << radius << "\n" << std::flush;

    Benchmark radius_bench;

    // Phase 1 (PARALLEL): Generate moves + apply SPR overlays
    // Each slot is independent — safe for concurrent access.
    using SPRStorageType = decltype(AddSPRStorageWithBackend<Backend>(sampled));
    std::vector<std::optional<SPRStorageType>> spr_slots(moves_per_radius);

    std::vector<size_t> work_items(moves_per_radius);
    std::iota(work_items.begin(), work_items.end(), 0);

    ParallelForEach(work_items, [&](size_t idx) {
      uint32_t thread_seed = seed ^ static_cast<uint32_t>(radius * 1000000 + idx);
      RandomMoveGenerator<SampledConstView> move_gen{sampled_const, thread_seed};

      auto move = move_gen.GenerateMove(radius);
      if (not move.has_value()) {
        return;
      }

      spr_slots[idx].emplace(AddSPRStorageWithBackend<Backend>(sampled));
      auto& spr = *spr_slots[idx];
      bool success =
          spr.View().InitHypotheticalTree(move->src, move->dst, move->lca);
      if (not success) {
        spr_slots[idx].reset();
      }
    });

    auto apply_ms = radius_bench.lapMs();

    // Phase 2 (PARALLEL): Create fragments and overlay MappedNodes.
    // Each fragment only reads its own SPR overlay and creates independent storage.
    // SetOverlay<MappedNodes> only writes to per-overlay replaced_node_storage_.
    using FragStorageType =
        decltype(std::declval<SPRStorageType&>().View().MakeFragment());
    std::vector<std::optional<FragStorageType>> frag_slots(moves_per_radius);

    ParallelForEach(work_items, [&](size_t idx) {
      if (not spr_slots[idx].has_value()) {
        return;
      }
      auto& spr = *spr_slots[idx];
      frag_slots[idx].emplace(spr.View().MakeFragment());
      for (auto node : frag_slots[idx]->View().GetNodes()) {
        auto spr_node = spr.View().Get(node.GetId());
        if (not spr_node.IsAppended()) {
          spr_node.template SetOverlay<MappedNodes>();
        }
      }
    });

    auto fragment_ms = radius_bench.lapMs();

    // Phase 3 (BATCH): Merge all fragments at once via AddDAGs.
    // AddDAGs internally parallelizes MergeCompactGenomes, MergeNodes, MergeEdges.
    // The reverse_map_ race in SetOriginalId is guarded by reverse_map_mtx_.
    using FragViewType = decltype(std::declval<FragStorageType&>().View());
    std::vector<FragViewType> frag_views;
    for (size_t idx = 0; idx < moves_per_radius; idx++) {
      if (frag_slots[idx].has_value()) {
        frag_views.push_back(frag_slots[idx]->View());
      }
    }
    if (not frag_views.empty()) {
      merge.AddDAGs(frag_views);
    }
    size_t accepted = frag_views.size();

    auto merge_ms = radius_bench.lapMs();

    std::cout << "  Accepted " << accepted << "/" << frag_views.size()
              << " generated (" << moves_per_radius
              << " attempted), apply=" << apply_ms
              << "ms frag=" << fragment_ms
              << "ms merge=" << merge_ms << "ms\n" << std::flush;

    results.push_back(RadiusResult{accepted, apply_ms, fragment_ms, merge_ms});
  }

  return results;
}

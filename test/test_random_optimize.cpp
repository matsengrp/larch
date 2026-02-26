#include "test_common.hpp"
#include "sample_dag.hpp"
#include "larch/spr/random_optimize.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"

// Test 1: Single random move — generate, apply via SPR, verify fragment is valid
static void test_single_random_move() {
  auto input_storage = make_sample_dag();
  auto dag = input_storage.View();

  // Set up merge with initial DAG
  Merge merge{dag.GetReferenceSequence()};
  merge.AddDAGs(std::vector{dag.Const()});
  merge.ComputeResultEdgeMutations();

  // Sample a tree
  SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
  auto sampled_storage = weight.MinWeightSampleTree({});
  auto sampled = sampled_storage.View();
  sampled.RecomputeCompactGenomes(true);
  sampled.SampleIdsFromCG();

  TestAssert(sampled.IsTree());

  // Generate a random move
  using SampledView = decltype(sampled.Const());
  RandomMoveGenerator<SampledView> move_gen{sampled.Const(), 42};

  auto move = move_gen.GenerateMove();
  TestAssert(move.has_value());

  // Apply the move via SPR with ParsimonyOnlyScoringBackend
  using SampledStorage = std::remove_reference_t<decltype(sampled_storage)>;
  using Backend = ParsimonyOnlyScoringBackend<SampledStorage>;
  auto spr = AddSPRStorageWithBackend<Backend>(sampled);

  bool success = spr.View().InitHypotheticalTree(move->src, move->dst, move->lca);
  TestAssert(success);

  // Extract fragment
  auto fragment = spr.View().MakeFragment();

  // Verify fragment has nodes and edges
  size_t node_count = 0;
  for ([[maybe_unused]] auto node : fragment.View().GetNodes()) {
    node_count++;
  }
  TestAssert(node_count > 0);

  size_t edge_count = 0;
  for ([[maybe_unused]] auto edge : fragment.View().GetEdges()) {
    edge_count++;
  }
  TestAssert(edge_count > 0);

  std::cout << "  Single move: src=" << move->src.value << " dst=" << move->dst.value
            << " lca=" << move->lca.value << " fragment_nodes=" << node_count
            << " fragment_edges=" << edge_count << "\n";
}

// Test 2: Run OptimizeDAGWithRandomMoves for a few iterations, verify DAG grows
static void test_random_optimize_grows_dag() {
  auto input_storage = make_sample_dag();
  auto dag = input_storage.View();

  Merge merge{dag.GetReferenceSequence()};
  merge.AddDAGs(std::vector{dag.Const()});

  size_t initial_nodes = merge.GetResult().GetNodesCount();
  size_t initial_edges = merge.GetResult().GetEdgesCount();

  std::cout << "  Initial: " << initial_nodes << " nodes, " << initial_edges
            << " edges\n";

  // Run optimization
  OptimizeDAGWithRandomMoves(merge, /*num_iterations=*/3, /*moves_per_iteration=*/10,
                             /*collapse_empty_fragment_edges=*/true,
                             /*seed=*/42);

  merge.ComputeResultEdgeMutations();
  size_t final_nodes = merge.GetResult().GetNodesCount();
  size_t final_edges = merge.GetResult().GetEdgesCount();

  std::cout << "  Final: " << final_nodes << " nodes, " << final_edges << " edges\n";

  // The DAG should have at least as many nodes/edges as before (merging adds structure)
  TestAssert(final_nodes >= initial_nodes);
  TestAssert(final_edges >= initial_edges);
}

// Test 3: Verify merged DAG is valid — edges have valid mutations, CGs are consistent
static void test_random_optimize_valid_dag() {
  auto input_storage = make_sample_dag();
  auto dag = input_storage.View();

  Merge merge{dag.GetReferenceSequence()};
  merge.AddDAGs(std::vector{dag.Const()});

  OptimizeDAGWithRandomMoves(merge, /*num_iterations=*/2, /*moves_per_iteration=*/5,
                             /*collapse_empty_fragment_edges=*/true,
                             /*seed=*/123);

  merge.ComputeResultEdgeMutations();
  auto result = merge.GetResult();

  // Verify internal->internal edges have consistent mutations.
  // (Leaf node CGs are empty in the merge result; edge mutations for leaves
  // are computed from the sample_id_to_cg_map, so we skip leaf children.)
  for (auto edge : result.GetEdges()) {
    auto parent = edge.GetParent();
    auto child = edge.GetChild();

    if (parent.IsUA() || child.IsLeaf()) {
      continue;
    }

    const auto& parent_cg = parent.GetCompactGenome();
    const auto& child_cg = child.GetCompactGenome();
    const auto& ref = result.GetReferenceSequence();

    for (auto& [pos, bases] : edge.GetEdgeMutations()) {
      auto parent_base = parent_cg.GetBase(pos, ref);
      auto child_base = child_cg.GetBase(pos, ref);
      TestAssert(parent_base == bases.first);
      TestAssert(child_base == bases.second);
    }
  }
}

// Test 4: Verify parsimony score is computable on result trees
static void test_random_optimize_parsimony_score() {
  auto input_storage = make_sample_dag();
  auto dag = input_storage.View();

  Merge merge{dag.GetReferenceSequence()};
  merge.AddDAGs(std::vector{dag.Const()});

  OptimizeDAGWithRandomMoves(merge, /*num_iterations=*/2, /*moves_per_iteration=*/5,
                             /*collapse_empty_fragment_edges=*/true,
                             /*seed=*/456);

  merge.ComputeResultEdgeMutations();

  // Sample a tree from the result and compute its parsimony score
  SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
  auto sampled = weight.MinWeightSampleTree({});
  TestAssert(sampled.View().IsTree());

  // Count edge mutations as a simple parsimony score proxy
  size_t total_mutations = 0;
  for (auto edge : sampled.View().GetEdges()) {
    if (not edge.GetParent().IsUA()) {
      for ([[maybe_unused]] auto& [pos, bases] : edge.GetEdgeMutations()) {
        total_mutations++;
      }
    }
  }
  std::cout << "  Sampled tree parsimony (mutation count): " << total_mutations << "\n";
  TestAssert(total_mutations > 0);
}

// Test 5: RandomMoveGenerator produces valid moves
static void test_random_move_generator() {
  auto input_storage = make_sample_dag();
  auto dag = input_storage.View();

  Merge merge{dag.GetReferenceSequence()};
  merge.AddDAGs(std::vector{dag.Const()});
  merge.ComputeResultEdgeMutations();

  SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
  auto sampled_storage = weight.MinWeightSampleTree({});
  auto sampled = sampled_storage.View();
  sampled.RecomputeCompactGenomes(true);
  sampled.SampleIdsFromCG();

  using SampledView = decltype(sampled.Const());
  RandomMoveGenerator<SampledView> move_gen{sampled.Const(), 99};

  size_t valid_moves = 0;
  for (size_t i = 0; i < 20; i++) {
    auto move = move_gen.GenerateMove();
    if (move.has_value()) {
      // Verify basic constraints
      TestAssert(move->src != move->dst);
      TestAssert(not sampled.Get(move->src).IsUA());
      TestAssert(not sampled.Get(move->dst).IsUA());
      valid_moves++;
    }
  }
  std::cout << "  Generated " << valid_moves << "/20 valid moves\n";
  TestAssert(valid_moves > 0);
}

[[maybe_unused]] static const auto test_added_rm1 =
    add_test({test_single_random_move, "random-moves: Single random move",
              {"random-moves"}});

[[maybe_unused]] static const auto test_added_rm2 =
    add_test({test_random_optimize_grows_dag, "random-moves: DAG grows with optimization",
              {"random-moves"}});

[[maybe_unused]] static const auto test_added_rm3 =
    add_test({test_random_optimize_valid_dag, "random-moves: Valid DAG after optimization",
              {"random-moves"}});

[[maybe_unused]] static const auto test_added_rm4 = add_test(
    {test_random_optimize_parsimony_score,
     "random-moves: Parsimony score computable", {"random-moves"}});

[[maybe_unused]] static const auto test_added_rm5 =
    add_test({test_random_move_generator, "random-moves: Move generator produces valid moves",
              {"random-moves"}});

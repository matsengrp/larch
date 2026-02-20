#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/spr/random_optimize.hpp"
#include "larch/benchmark.hpp"

#include <sys/resource.h>

static MADAGStorage<> Load(std::string_view input_dag_path,
                           std::string_view refseq_path) {
  std::string reference_sequence = LoadReferenceSequence(refseq_path);
  MADAGStorage<> input_dag_storage =
      LoadTreeFromProtobuf(input_dag_path, reference_sequence);
  input_dag_storage.View().RecomputeCompactGenomes(true);
  return input_dag_storage;
}

static MADAGStorage<> Load(std::string_view input_dag_path) {
  MADAGStorage<> input_dag_storage = LoadDAGFromProtobuf(input_dag_path);
  input_dag_storage.View().RecomputeCompactGenomes(true);
  input_dag_storage.View().SampleIdsFromCG(true);
  return input_dag_storage;
}

static void test_spr_random(const MADAGStorage<>& input_dag_storage, size_t count,
                             size_t moves_per_radius = 100) {
  MADAG input_dag = input_dag_storage.View();
  Merge merge{input_dag.GetReferenceSequence()};
  merge.AddDAGs(std::vector{input_dag});

  std::mt19937 seed_gen{42};

  for (size_t i = 0; i < count; ++i) {
    Benchmark phase_bench;

    // 1. Compute edge mutations
    std::cout << "  Computing edge mutations..." << std::flush;
    merge.ComputeResultEdgeMutations();
    auto compute_edge_ms = phase_bench.lapMs();
    std::cout << " done (" << compute_edge_ms << "ms)\n" << std::flush;

    // 2. Sample a min-weight tree
    std::cout << "  Sampling tree..." << std::flush;
    SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
    auto sampled_storage = weight.MinWeightSampleTree({});
    auto sampled = sampled_storage.View();
    sampled.RecomputeCompactGenomes(true);
    sampled.SampleIdsFromCG();
    Assert(sampled.IsTree());
    auto sampling_ms = phase_bench.lapMs();
    std::cout << " done (" << sampling_ms << "ms)\n" << std::flush;

    // 3. Parallel optimization with radius iteration
    uint32_t iter_seed = seed_gen();
    auto radius_results = OptimizeDAGParallelRadius(merge, sampled_storage,
                                                     moves_per_radius, iter_seed);
    auto optimize_ms = phase_bench.lapMs();

    // 4. Merge the sampled tree itself
    merge.AddDAG(sampled);
    auto merge_tree_ms = phase_bench.lapMs();

    // 5. Report parsimony score
    merge.ComputeResultEdgeMutations();
    SubtreeWeight<ParsimonyScore, MergeDAG> post_weight{merge.GetResult()};
    auto best = post_weight.MinWeightSampleTree({});
    best.View().RecomputeCompactGenomes(true);
    size_t parsimony = 0;
    for (auto edge : best.View().GetEdges()) {
      if (not edge.GetParent().IsUA()) {
        for ([[maybe_unused]] auto& [pos, bases] : edge.GetEdgeMutations()) {
          parsimony++;
        }
      }
    }
    auto score_report_ms = phase_bench.lapMs();

    // Summarize radius results
    long apply_total_ms = 0;
    long fragment_total_ms = 0;
    long merge_total_ms = 0;
    size_t total_accepted = 0;
    for (auto& r : radius_results) {
      apply_total_ms += r.apply_ms;
      fragment_total_ms += r.fragment_ms;
      merge_total_ms += r.merge_ms;
      total_accepted += r.accepted_moves;
    }

    auto serial_outside_optimize =
        compute_edge_ms + sampling_ms + merge_tree_ms + score_report_ms;
    auto total_ms =
        compute_edge_ms + sampling_ms + optimize_ms + merge_tree_ms + score_report_ms;

    std::cout << "=== Iteration " << (i + 1) << " timing (ms) ===\n"
              << "  [SERIAL]   ComputeEdgeMutations: " << compute_edge_ms << "\n"
              << "  [SERIAL]   Sampling:             " << sampling_ms << "\n"
              << "  [MIXED]    Optimize (" << total_accepted << " moves, "
              << radius_results.size() << " radii): " << optimize_ms << "\n"
              << "    [PARALLEL] Apply total:         " << apply_total_ms << " ms\n"
              << "    [PARALLEL] Fragment total:      " << fragment_total_ms << " ms\n"
              << "    [BATCH]    Merge total:         " << merge_total_ms << " ms\n"
              << "  [SERIAL]   MergeSampledTree:     " << merge_tree_ms << "\n"
              << "  [SERIAL]   ScoreReport:          " << score_report_ms << "\n"
              << "  Serial (outside optimize): " << serial_outside_optimize << " ms\n"
              << "  Total: " << total_ms << " ms\n";

    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "Parsimony score after iteration " << (i + 1) << ": " << parsimony
              << " (DAG nodes=" << merge.GetResult().GetNodesCount()
              << " edges=" << merge.GetResult().GetEdgesCount() << ")"
              << " peak_mem=" << (usage.ru_maxrss / 1024) << "MB\n" << std::flush;
  }
}

[[maybe_unused]] static const auto test_added_srm0 =
    add_test({[] {
                auto input = Load("data/test_5_trees/tree_0.pb.gz");
                Benchmark bench;
                bench.start();
                test_spr_random(input, 3, 100);
                bench.stop();
                std::cout << "Total time: " << bench.durationMs() << " ms\n";
              },
              "SPR-random: test_5_trees"});

[[maybe_unused]] static const auto test_added_srm1 =
    add_test({[] {
                auto input =
                    Load("data/seedtree/seedtree.pb.gz", "data/seedtree/refseq.txt.gz");
                Benchmark bench;
                bench.start();
                test_spr_random(input, 1, 10);
                bench.stop();
                std::cout << "Total time: " << bench.durationMs() << " ms\n";
              },
              "SPR-random: seedtree-quick"});

[[maybe_unused]] static const auto test_added_srm2 =
    add_test({[] {
                auto input =
                    Load("data/seedtree/seedtree.pb.gz", "data/seedtree/refseq.txt.gz");
                Benchmark bench;
                bench.start();
                test_spr_random(input, 3, 1000);
                bench.stop();
                std::cout << "Total time: " << bench.durationMs() << " ms\n";
              },
              "SPR-random: seedtree",
              {"slow"}});

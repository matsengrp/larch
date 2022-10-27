#include "mutation_annotated_dag.hpp"
#include "subtree_weight.hpp"
#include "parsimony_score.hpp"

#include <iostream>
#include <string_view>

#include "test_common.hpp"
#include "dag_loader.hpp"
#include "tree_count.hpp"
#include "merge.hpp"
#include "benchmark.hpp"

#if __has_include(<valgrind/callgrind.h>)
#include <valgrind/callgrind.h>
#endif

static void test_sample_tree(MADAG& dag) {
  if (dag.GetEdgeMutations().empty()) {
    dag.RecomputeEdgeMutations();
  }

  SubtreeWeight<ParsimonyScore> weight(dag);

  MADAG result = weight.SampleTree({}).first;
  assert_true(result.GetDAG().IsTree(), "Tree");

  SubtreeWeight<TreeCount> tree_count{dag};
  MADAG result2 = tree_count.UniformSampleTree({}).first;
  assert_true(result2.GetDAG().IsTree(), "Tree");
}

static void test_sample_tree(std::string_view path) {
  MADAG dag = LoadDAGFromProtobuf(path);
  test_sample_tree(dag);
}

static void bench_sampling(std::string_view path, std::string_view refseq_path) {
  MADAG dag = LoadTreeFromProtobuf(path, LoadReferenceSequence(refseq_path));
#if defined(CALLGRIND_START_INSTRUMENTATION)
  CALLGRIND_START_INSTRUMENTATION;
#endif
  Benchmark bench;
  bench.start();
for (int i = 0; i < 100; ++i) {
  Merge merge{dag.GetReferenceSequence()};
  merge.AddDAGs({dag});
  merge.ComputeResultEdgeMutations();
  SubtreeWeight<ParsimonyScore> weight{merge.GetResult()};
  std::ignore = weight.SampleTree({});
}
  bench.stop();
#if defined(CALLGRIND_START_INSTRUMENTATION)
  CALLGRIND_STOP_INSTRUMENTATION;
#endif
  std::cout << " " << bench.durationMs() << " ms ";
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_sample_tree("data/testcase/full_dag.pb.gz"); },
              "Sample tree: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_sample_tree("data/testcase1/full_dag.pb.gz"); },
              "Sample tree: testcase1"});

[[maybe_unused]] static const auto test_added2 =
    add_test({[] {
                bench_sampling("data/AY.103/AY.103_start_tree_no_ancestral.pb.gz",
                               "data/AY.103/ref_seq_noancestral.txt.gz");
              },
              "Bench sample tree: AY.103"});

[[maybe_unused]] static const auto test_added3 =
    add_test({[] {
                bench_sampling("data/20B/20B_start_tree_no_ancestral.pb.gz",
                               "data/20B/ref_seq_noancestral.txt.gz");
              },
              "Bench sample tree: 20B"});

[[maybe_unused]] static const auto test_added4 =
    add_test({[] {
                bench_sampling("data/B.1.1.529/B.1.1.529_start_tree_no_ancestral.pb.gz",
                               "data/B.1.1.529/ref_seq_noancestral.txt.gz");
              },
              "Bench sample tree: B.1.1.529"});
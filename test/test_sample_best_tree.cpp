#include "larch/madag/mutation_annotated_dag.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/subtree/parsimony_score_binary.hpp"

#include <iostream>
#include <string_view>

#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/subtree/tree_count.hpp"
#include "larch/merge/merge.hpp"
#include "larch/benchmark.hpp"

#if __has_include(<valgrind/callgrind.h>)
#include <valgrind/callgrind.h>
#endif

static void test_sample_tree(MADAG dag) {
  SubtreeWeight<ParsimonyScore, MADAG> weight(dag);

  auto result = weight.SampleTree({});
  TestAssert(result.View().IsTree());

  SubtreeWeight<TreeCount, MADAG> tree_count{dag};
  auto result2 = tree_count.UniformSampleTree({});
  TestAssert(result2.View().IsTree());

  SubtreeWeight<BinaryParsimonyScore, MADAG> binary_weight{dag};
  auto result3 = binary_weight.SampleTree({});
  TestAssert(result3.View().IsTree());

  Merge merge{dag.GetReferenceSequence()};
  merge.AddDAGs(std::vector{dag});
  merge.GetResult().GetRoot().Validate(true, true);
  merge.ComputeResultEdgeMutations();

  SubtreeWeight<BinaryParsimonyScore, MergeDAG> binary_merge_weight{merge.GetResult()};
  auto result4 = binary_merge_weight.SampleTree({});
  TestAssert(result4.View().IsTree());
}

static void test_sample_tree(std::string_view path) {
  MADAGStorage dag = LoadDAGFromProtobuf(path);
  dag.View().RecomputeCompactGenomes(true);
  dag.View().SampleIdsFromCG(true);
  test_sample_tree(dag.View());
}

[[maybe_unused]] static void bench_sampling(std::string_view path,
                                            std::string_view refseq_path) {
  MADAGStorage dag = LoadTreeFromProtobuf(path, LoadReferenceSequence(refseq_path));
  dag.View().RecomputeCompactGenomes(true);
#if defined(CALLGRIND_START_INSTRUMENTATION)
  CALLGRIND_START_INSTRUMENTATION;
#endif
  Benchmark bench;
  bench.start();
  Merge merge{dag.View().GetReferenceSequence()};
  merge.AddDAGs(std::vector{dag.View()});
  merge.ComputeResultEdgeMutations();
  SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
  std::ignore = weight.MinWeightSampleTree({});
  bench.stop();
#if defined(CALLGRIND_START_INSTRUMENTATION)
  CALLGRIND_STOP_INSTRUMENTATION;
#endif
  std::cout << " " << bench.durationMs() << " ms ";
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_sample_tree("data/testcase/full_dag.pb.gz"); },
              "Sample min weight tree: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_sample_tree("data/testcase1/full_dag.pb.gz"); },
              "Sample min weight tree: testcase1"});

// [[maybe_unused]] static const auto test_added2 =
//     add_test({[] {
//                 bench_sampling("data/AY.103/AY.103_start_tree_no_ancestral.pb.gz",
//                                "data/AY.103/ref_seq_noancestral.txt.gz");
//               },
//               "Bench sample min weight tree: AY.103"});

// [[maybe_unused]] static const auto test_added3 =
//     add_test({[] {
//                 bench_sampling("data/20B/20B_start_tree_no_ancestral.pb.gz",
//                                "data/20B/ref_seq_noancestral.txt.gz");
//               },
//               "Bench sample min weight tree: 20B"});

// [[maybe_unused]] static const auto test_added4 =
//     add_test({[] {
//                 bench_sampling("data/B.1.1.529/B.1.1.529_start_tree_no_ancestral.pb.gz",
//                                "data/B.1.1.529/ref_seq_noancestral.txt.gz");
//               },
//               "Bench sample min weight tree: B.1.1.529"});

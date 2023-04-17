#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score_binary.hpp"
#include "larch/subtree/tree_count.hpp"

#include <iostream>
#include <string_view>
#include <vector>

#include "test_common.hpp"

#include "larch/dag_loader.hpp"

static void test_dag_trimming(MADAG dag) {
  SubtreeWeight<BinaryParsimonyScore, MADAG> weight{dag};
  SubtreeWeight<TreeCount, MADAG> treecount{dag};
  auto before_trim_count = treecount.ComputeWeightBelow(dag.GetRoot(), {});
  auto before_trim_optimal_count = weight.MinWeightCount(dag.GetRoot(), {});

  MADAGStorage trimmed = weight.TrimToMinWeight({});
  MADAG trimmed_view = trimmed.View();
  SubtreeWeight<BinaryParsimonyScore, MADAG> trimmed_weight{trimmed_view};
  auto after_trim_optimal_count = trimmed_weight.MinWeightCount(trimmed_view.GetRoot(), {});
  SubtreeWeight<TreeCount, MADAG> trimtreecount{trimmed_view};
  auto after_trim_count = trimtreecount.ComputeWeightBelow(trimmed_view.GetRoot(), {});

  std::cout << "Before trim " << before_trim_optimal_count << "/" << before_trim_count << "optimal\n" << std::flush;
  std::cout << "After trim " << after_trim_optimal_count << "/" << after_trim_count << "optimal\n" << std::flush;

  assert_equal(after_trim_optimal_count, before_trim_optimal_count, "Before and after trim optimal counts");
  assert_equal(after_trim_count, before_trim_optimal_count, "Optimal trees before and all trees after trim");
}

static void test_dag_trimming(std::string_view path) {
  MADAGStorage dag = LoadDAGFromProtobuf(path);
  test_dag_trimming(dag.View());
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_dag_trimming("data/testcase/full_dag.pb.gz"); },
              "DAG trimming: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_dag_trimming("data/testcase1/full_dag.pb.gz"); },
              "DAG trimming: testcase1"});

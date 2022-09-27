#include "subtree_weight.hpp"
#include "parsimony_score.hpp"

#include <iostream>
#include <string_view>
#include <vector>

#include "test_common.hpp"

#include "dag_loader.hpp"

static void test_dag_trimming(MADAG& dag, size_t expected_edges) {
  if (dag.GetEdgeMutations().empty()) {
    dag.RecomputeEdgeMutations();
  }

  SubtreeWeight<ParsimonyScore> weight{dag};

  MADAG trimmed = weight.TrimToMinWeight({});

  assert_equal(trimmed.GetDAG().GetEdgesCount(), expected_edges, "Edges count");
}

static void test_dag_trimming(std::string_view path, size_t expected_edges) {
  MADAG dag = LoadDAGFromProtobuf(path);
  test_dag_trimming(dag, expected_edges);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_dag_trimming("data/testcase/full_dag.pb.gz", 57); },
              "DAG trimming: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_dag_trimming("data/testcase1/full_dag.pb.gz", 59); },
              "DAG trimming: testcase1"});

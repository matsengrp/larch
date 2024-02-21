#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score_binary.hpp"

#include <iostream>
#include <string_view>
#include <vector>

#include "test_common.hpp"

#include "larch/dag_loader.hpp"

static void test_dag_trimming(MADAG dag, size_t expected_edges) {
  SubtreeWeight<BinaryParsimonyScore, MADAG> weight{dag};

  MADAGStorage trimmed = weight.TrimToMinWeight({});

  TestAssert(trimmed.View().GetEdgesCount() == expected_edges);
}

static void test_dag_trimming(std::string_view path, size_t expected_edges) {
  MADAGStorage dag = LoadDAGFromProtobuf(path);
  test_dag_trimming(dag.View(), expected_edges);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_dag_trimming("data/testcase/full_dag.pb.gz", 530); },
              "DAG trimming: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_dag_trimming("data/testcase1/full_dag.pb.gz", 59); },
              "DAG trimming: testcase1"});

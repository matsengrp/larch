#include "subtree_weight.hpp"
#include "parsimony_score.hpp"

#include <iostream>
#include <string_view>
#include <vector>

#include "test_common.hpp"

#include "dag_loader.hpp"

static void test_dag_trimming(MADAG& dag, size_t expected_edges) {
  if (dag.GetEdgeMutations().empty()) {
    dag.GetEdgeMutations() = dag.ComputeEdgeMutations(dag.GetReferenceSequence());
  }

  SubtreeWeight<ParsimonyScore> weight{dag};

  MADAG trimmed = weight.TrimToMinWeight({});

  size_t edges_count = 0;
  for (Edge edge : trimmed.GetDAG().GetEdges()) {
    if (not(edge.GetId() == EdgeId{})) {
      ++edges_count;  // Manually counting valid edges until reindexing is ready
    }
  }

  assert_equal(edges_count, expected_edges, "Edgs count");
}

static void test_dag_trimming(std::string_view path, size_t expected_edges) {
  MADAG dag = LoadDAGFromProtobuf(path);
  test_dag_trimming(dag, expected_edges);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_dag_trimming("data/testcase/full_dag.pb.gz", 2308); },
              "DAG trimming: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_dag_trimming("data/testcase1/full_dag.pb.gz", 203); },
              "DAG trimming: testcase1"});

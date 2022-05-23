#include "merge.hpp"

#include <iostream>
#include <fstream>

#include "test_common.hpp"
#include "history_dag_loader.hpp"
#include "benchmark.hpp"

static void test_protobuf(const std::string& correct_path,
                          const std::vector<std::string>& paths) {
  std::vector<std::vector<Mutations>> mutations;
  std::vector<HistoryDAG> trees;
  for (auto& path : paths) {
    std::vector<Mutations> tree_mutations;
    std::string ref_seq;
    trees.emplace_back(LoadHistoryDAGFromProtobufGZ(path, ref_seq, tree_mutations));
    mutations.emplace_back(std::move(tree_mutations));
  }
  std::string reference_sequence;
  HistoryDAG correct_result =
      LoadHistoryDAGFromJsonGZ(correct_path, reference_sequence);

  Merge merge(reference_sequence, trees, mutations);
  merge.Run();

  assert_equal(correct_result.GetNodes().size(), merge.GetResult().GetNodes().size(),
               "Nodes count");

  assert_equal(correct_result.GetEdges().size(), merge.GetResult().GetEdges().size(),
               "Edges count");
}

static void test_five_trees() {
  test_protobuf("data/test_5_trees/full_dag.json.gz",
                {
                    "data/test_5_trees/tree_0.pb.gz",
                    "data/test_5_trees/tree_1.pb.gz",
                    "data/test_5_trees/tree_2.pb.gz",
                    "data/test_5_trees/tree_3.pb.gz",
                    "data/test_5_trees/tree_4.pb.gz",
                });
}

static void test_case_2() {
  test_protobuf("data/testcase2/full_dag.json.gz", {
                                                       "data/testcase2/tree_0.pb.gz",
                                                       "data/testcase2/tree_1.pb.gz",
                                                       "data/testcase2/tree_2.pb.gz",
                                                       "data/testcase2/tree_3.pb.gz",
                                                       "data/testcase2/tree_4.pb.gz",
                                                   });
}

static void test_case_ref() {
  test_protobuf("data/testcaseref/full_dag.json.gz",
                {
                    "data/testcaseref/tree_0.pb.gz",
                    "data/testcaseref/tree_1.pb.gz",
                    "data/testcaseref/tree_2.pb.gz",
                    "data/testcaseref/tree_3.pb.gz",
                    "data/testcaseref/tree_4.pb.gz",
                });
}

[[maybe_unused]] static const auto test0_added = add_test({test_case_2, "Test case 2"});

[[maybe_unused]] static const auto test1_added =
    add_test({test_five_trees, "Merge 5 trees"});

[[maybe_unused]] static const auto test2_added =
    add_test({test_case_ref, "Tree with different ref"});

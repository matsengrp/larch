#include "merge.hpp"

#include <iostream>
#include <fstream>
#include <experimental/filesystem>

#include "test_common.hpp"
#include "dag_loader.hpp"
#include "benchmark.hpp"

static void test_protobuf(const std::string& correct_path,
                          const std::vector<std::string>& paths) {
  std::vector<MADAG> trees;
  for (auto& path : paths) {
    trees.emplace_back(LoadDAGFromProtobuf(path));
  }

  MADAG correct_result = LoadDAGFromJson(correct_path);

  Merge merge(correct_result.reference_sequence);
  std::vector<std::reference_wrapper<MADAG>> tree_refs{trees.begin(), trees.end()};
  merge.AddTrees(tree_refs);

  assert_equal(correct_result.dag.GetNodesCount(), merge.GetResult().GetNodesCount(),
               "Nodes count");

  assert_equal(correct_result.dag.GetEdgesCount(), merge.GetResult().GetEdgesCount(),
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

static void test_case_20d() {
  std::vector<std::string> paths;
  const std::experimental::filesystem::path dir{"data/20D_from_fasta"};
  for (auto& file : std::experimental::filesystem::directory_iterator{dir}) {
    std::string name = file.path().filename().string();
    if (name.find("1final-tree-") == 0) {
      paths.push_back(file.path().string());
    }
  }

  std::vector<MADAG> trees;

  MADAG correct_result = LoadDAGFromJson("data/20D_from_fasta/20D_full_dag.json.gz");

  trees.resize(paths.size());
  std::vector<std::pair<size_t, std::string_view>> paths_idx;
  for (size_t i = 0; i < paths.size(); ++i) {
    paths_idx.push_back({i, paths.at(i)});
  }
  tbb::parallel_for_each(paths_idx.begin(), paths_idx.end(), [&](auto path_idx) {
    trees.at(path_idx.first) = LoadTreeFromProtobuf(path_idx.second);
  });

  Benchmark merge_time;
  Merge merge(correct_result.reference_sequence);
  std::vector<std::reference_wrapper<MADAG>> tree_refs{trees.begin(), trees.end()};
  merge_time.start();
  merge.AddTrees(tree_refs, false);
  merge_time.stop();
  std::cout << " DAGs merged in " << merge_time.durationMs() << " ms. ";

  assert_equal(correct_result.dag.GetNodesCount(), merge.GetResult().GetNodesCount(),
               "Nodes count");

  assert_equal(correct_result.dag.GetEdgesCount(), merge.GetResult().GetEdgesCount(),
               "Edges count");
}

static void test_add_trees() {
  std::string_view correct_path = "data/test_5_trees/full_dag.json.gz";
  std::vector<std::string> paths1 = {"data/test_5_trees/tree_0.pb.gz",
                                     "data/test_5_trees/tree_1.pb.gz"};
  std::vector<std::string> paths2 = {"data/test_5_trees/tree_2.pb.gz",
                                     "data/test_5_trees/tree_3.pb.gz",
                                     "data/test_5_trees/tree_4.pb.gz"};

  std::vector<MADAG> trees1, trees2;
  for (auto& path : paths1) {
    trees1.emplace_back(LoadDAGFromProtobuf(path));
  }
  for (auto& path : paths2) {
    trees2.emplace_back(LoadDAGFromProtobuf(path));
  }

  MADAG correct_result = LoadDAGFromJson(correct_path);

  Merge merge(correct_result.reference_sequence);
  std::vector<std::reference_wrapper<MADAG>> tree_refs1{trees1.begin(), trees1.end()};
  merge.AddTrees(tree_refs1);
  std::vector<std::reference_wrapper<MADAG>> tree_refs2{trees2.begin(), trees2.end()};
  merge.AddTrees(tree_refs2);

  assert_equal(correct_result.dag.GetNodesCount(), merge.GetResult().GetNodesCount(),
               "Nodes count");

  assert_equal(correct_result.dag.GetEdgesCount(), merge.GetResult().GetEdgesCount(),
               "Edges count");
}

[[maybe_unused]] static const auto test0_added =
    add_test({test_case_2, "Merge: Test case 2"});

[[maybe_unused]] static const auto test1_added =
    add_test({test_five_trees, "Merge: 5 trees"});

[[maybe_unused]] static const auto test2_added =
    add_test({test_case_ref, "Merge: Tree with different ref"});

[[maybe_unused]] static const auto test3_added =
    add_test({test_case_20d, "Merge: 800 trees"});

[[maybe_unused]] static const auto test4_added =
    add_test({test_add_trees, "Merge: Add trees"});

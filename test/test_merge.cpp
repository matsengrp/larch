#include "merge.hpp"

#include <iostream>
#include <fstream>
#include <experimental/filesystem>

#include "test_common.hpp"
#include "history_dag_loader.hpp"
#include "benchmark.hpp"

#include <valgrind/callgrind.h>  //XXX

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

static void test_case_20d() {
  std::vector<std::string> paths;
  const std::experimental::filesystem::path dir{"data/20D_from_fasta"};
  for (auto& file : std::experimental::filesystem::directory_iterator{dir}) {
    if (file.path().filename().string() != "20D_full_dag.json.gz") {
      paths.push_back(file.path().string());
    }
  }

  std::vector<std::vector<Mutations>> mutations;
  std::vector<HistoryDAG> trees;
  std::string reference_sequence;

  HistoryDAG correct_result = LoadHistoryDAGFromJsonGZ(
      "data/20D_from_fasta/20D_full_dag.json.gz", reference_sequence);

  trees.resize(paths.size());
  mutations.resize(paths.size());
  std::vector<std::pair<size_t, std::string_view>> paths_idx;
  for (size_t i = 0; i < paths.size(); ++i) {
    paths_idx.push_back({i, paths.at(i)});
  }
  std::cout << "Loading trees ";
  tbb::parallel_for_each(paths_idx.begin(), paths_idx.end(), [&](auto path_idx) {
    std::vector<Mutations> tree_mutations;
    std::cout << "." << std::flush;
    trees.at(path_idx.first) = LoadTreeFromProtobufGZ(path_idx.second, tree_mutations);
    mutations.at(path_idx.first) = std::move(tree_mutations);
  });
  std::cout << " done."
            << "\n";

  Benchmark merge_time;
  Merge merge(reference_sequence, trees, mutations, true);
  merge_time.start();
  CALLGRIND_START_INSTRUMENTATION;
  merge.Run();
  merge_time.stop();
  std::cout << "\nDAGs merged in " << merge_time.durationMs() << " ms\n";

  std::cout << "Nodes: " << merge.GetResult().GetNodes().size() << "\n";
  std::cout << "Edges: " << merge.GetResult().GetEdges().size() << "\n";

  assert_equal(correct_result.GetNodes().size(), merge.GetResult().GetNodes().size(),
               "Nodes count");

  assert_equal(correct_result.GetEdges().size(), merge.GetResult().GetEdges().size(),
               "Edges count");
}

[[maybe_unused]] static const auto test0_added = add_test({test_case_2, "Test case 2"});

[[maybe_unused]] static const auto test1_added =
    add_test({test_five_trees, "Merge 5 trees"});

[[maybe_unused]] static const auto test2_added =
    add_test({test_case_ref, "Tree with different ref"});

[[maybe_unused]] static const auto test3_added = add_test({test_case_20d, "800 trees"});

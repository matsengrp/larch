#include "larch/merge/merge.hpp"

#include <iostream>
#include <fstream>
#include <experimental/filesystem>
#include <vector>

#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "benchmark.hpp"

static void test_protobuf(const std::string& correct_path,
                          const std::vector<std::string>& paths) {
  std::vector<MADAGStorage> trees;
  trees.reserve(paths.size());
  std::vector<MADAG> tree_views;
  for (auto& path : paths) {
    trees.emplace_back(LoadDAGFromProtobuf(path));
    MutableMADAG view = trees.back().View();
    view.RecomputeCompactGenomes(true);
    view.SampleIdsFromCG();
    for (auto leaf : view.GetLeafs()) {
      Assert(leaf.HaveSampleId());
    }
    tree_views.push_back(view);
  }

  MADAGStorage correct_result = LoadDAGFromJson(correct_path);
  correct_result.View().RecomputeEdgeMutations();

  Merge merge(correct_result.View().GetReferenceSequence());
  merge.AddDAGs(tree_views);
  merge.GetResult().GetRoot().Validate(true, true);

  assert_equal(correct_result.View().GetNodesCount(), merge.GetResult().GetNodesCount(),
               "Nodes count");

  assert_equal(correct_result.View().GetEdgesCount(), merge.GetResult().GetEdgesCount(),
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

  std::vector<MADAGStorage> trees;

  MADAGStorage correct_result = LoadDAGFromJson("data/20D_from_fasta/full_dag.json.gz");
  correct_result.View().RecomputeEdgeMutations();

  trees.reserve(paths.size());
  std::vector<std::pair<size_t, std::string_view>> paths_idx;
  for (size_t i = 0; i < paths.size(); ++i) {
    trees.push_back(MADAGStorage{{}});
    paths_idx.push_back({i, paths.at(i)});
  }
  tbb::parallel_for_each(paths_idx.begin(), paths_idx.end(), [&](auto path_idx) {
    trees.at(path_idx.first) = LoadTreeFromProtobuf(
        path_idx.second, correct_result.View().GetReferenceSequence());
    for (auto node : trees.at(path_idx.first).View().GetNodes()) {
      if (node.IsLeaf()) {
        Assert(node.GetSampleId().has_value());
        Assert(not node.GetSampleId().value().empty());
      }
    }
    trees.at(path_idx.first).View().RecomputeCompactGenomes(true);
  });

  std::vector<MADAG> tree_views;
  for (auto& i : trees) {
    tree_views.emplace_back(i);
  }

  Benchmark merge_time;
  Merge merge(correct_result.View().GetReferenceSequence());
  merge_time.start();
  merge.AddDAGs(tree_views);
  merge_time.stop();
  merge.GetResult().GetRoot().Validate(true, true);
  std::cout << " DAGs merged in " << merge_time.durationMs() << " ms. ";

  assert_equal(correct_result.View().GetNodesCount(), merge.GetResult().GetNodesCount(),
               "Nodes count");

  assert_equal(correct_result.View().GetEdgesCount(), merge.GetResult().GetEdgesCount(),
               "Edges count");
}

static void test_add_trees() {
  std::string_view correct_path = "data/test_5_trees/full_dag.json.gz";
  std::vector<std::string> paths1 = {"data/test_5_trees/tree_0.pb.gz",
                                     "data/test_5_trees/tree_1.pb.gz"};
  std::vector<std::string> paths2 = {"data/test_5_trees/tree_2.pb.gz",
                                     "data/test_5_trees/tree_3.pb.gz",
                                     "data/test_5_trees/tree_4.pb.gz"};

  std::vector<MADAGStorage> trees1, trees2;
  for (auto& path : paths1) {
    trees1.push_back(LoadDAGFromProtobuf(path));
    trees1.back().View().RecomputeCompactGenomes(true);
    trees1.back().View().SampleIdsFromCG();
  }
  for (auto& path : paths2) {
    trees2.push_back(LoadDAGFromProtobuf(path));
    trees2.back().View().RecomputeCompactGenomes(true);
    trees2.back().View().SampleIdsFromCG();
  }

  MADAGStorage correct_result = LoadDAGFromJson(correct_path);
  correct_result.View().RecomputeEdgeMutations();

  Merge merge(correct_result.View().GetReferenceSequence());
  std::vector<MADAG> tree_views1;
  for (auto& i : trees1) {
    tree_views1.emplace_back(i);
  }
  merge.AddDAGs(tree_views1);
  merge.GetResult().GetRoot().Validate(true, true);
  std::vector<MADAG> tree_views2;
  for (auto& i : trees2) {
    tree_views2.emplace_back(i);
  }
  merge.AddDAGs(tree_views2);
  merge.GetResult().GetRoot().Validate(true, true);

  assert_equal(correct_result.View().GetNodesCount(), merge.GetResult().GetNodesCount(),
               "Nodes count");

  assert_equal(correct_result.View().GetEdgesCount(), merge.GetResult().GetEdgesCount(),
               "Edges count");
}

static void test_subtree() {
  std::string_view correct_path = "data/test_5_trees/full_dag.json.gz";
  std::vector<std::string> paths = {
      "data/test_5_trees/tree_0.pb.gz", "data/test_5_trees/tree_1.pb.gz",
      "data/test_5_trees/tree_2.pb.gz", "data/test_5_trees/tree_3.pb.gz",
      "data/test_5_trees/tree_4.pb.gz"};

  MADAGStorage correct_result = LoadDAGFromJson(correct_path);
  correct_result.View().RecomputeEdgeMutations();
  Merge merge(correct_result.View().GetReferenceSequence());

  std::vector<MADAGStorage> trees;
  for (auto& path : paths) {
    trees.push_back(LoadDAGFromProtobuf(path));
    trees.back().View().RecomputeCompactGenomes(true);
    trees.back().View().SampleIdsFromCG();
  }

  for (auto& tree : trees) {
    Fragment frag{tree.View(),
                  tree.View().GetNodes() | Transform::GetId() | ranges::to_vector,
                  tree.View().GetEdges() | Transform::GetId() | ranges::to_vector};
    merge.AddDAG(frag);
    merge.GetResult().GetRoot().Validate(true, true);
  }

  assert_equal(correct_result.View().GetNodesCount(), merge.GetResult().GetNodesCount(),
               "Nodes count");

  assert_equal(correct_result.View().GetEdgesCount(), merge.GetResult().GetEdgesCount(),
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

[[maybe_unused]] static const auto test5_added =
    add_test({test_subtree, "Merge: Subtree"});

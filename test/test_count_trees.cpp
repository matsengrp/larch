#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/tree_count.hpp"

#include <iostream>
#include <string_view>
#include <string>

#include "test_common.hpp"

#include "larch/dag_loader.hpp"

static void test_tree_count(MADAG dag, TreeCount::Weight expected_score) {
  SubtreeWeight<MADAG, TreeCount> treecount(dag);

  TreeCount::Weight score = treecount.ComputeWeightBelow(dag.GetRoot(), {});

  assert_equal(score, expected_score, "Tree count doesn't match truth");
}

static void test_tree_count(std::string_view path, TreeCount::Weight expected_score) {
  MADAGStorage dag = LoadDAGFromProtobuf(path);

  test_tree_count(dag.View(), expected_score);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_tree_count("data/testcase/full_dag.pb.gz", 818); },
              "Subtree weight: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_tree_count("data/testcase1/full_dag.pb.gz", 7); },
              "Subtree weight: testcase1"});

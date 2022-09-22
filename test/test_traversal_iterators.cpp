#include "test_common.hpp"
#include "dag.hpp"
#include "synthetic_dags.hpp"
#include "dag_loader.hpp"

static const std::vector<size_t> correct_preorder = {
    0,  1, 2, 4,  7,  12, 13, 8,  14, 15, 22, 27, 28, 29, 30, 23,
    16, 5, 9, 17, 18, 24, 25, 19, 11, 20, 26, 21, 6,  10, 3};

static const std::vector<size_t> correct_postorder = {
    12, 13, 7,  14, 27, 29, 30, 28, 22, 23, 15, 16, 8, 4, 17, 24,
    25, 18, 19, 9,  26, 20, 21, 11, 5,  10, 6,  2,  3, 1, 0};

[[maybe_unused]] static void test_traversal_iterators() {
  MADAG dag = MakeSyntheticDAG();
  auto correct_pre = correct_preorder.begin();
  for (Node node : dag.GetDAG().TraversePreOrder()) {
    assert_equal(node.GetId().value, *correct_pre++, "Preorder");
  }
  auto correct_post = correct_postorder.begin();
  for (Node node : dag.GetDAG().TraversePostOrder()) {
    assert_equal(node.GetId().value, *correct_post++, "Postorder");
  }
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_traversal_iterators(); }, "Traversal iterators"});

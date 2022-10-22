#include "test_common.hpp"

#include "dag_loader.hpp"

static void AssertDAGsEqual(const MADAG& lhs, const MADAG& rhs) {
  assert_equal(lhs.GetDAG().GetNodesCount(), rhs.GetDAG().GetNodesCount(),
               "Nodes count");
  assert_equal(lhs.GetDAG().GetEdgesCount(), rhs.GetDAG().GetEdgesCount(),
               "Edges count");
  auto lhs_order = lhs.GetDAG().TraversePostOrder();
  auto rhs_order = rhs.GetDAG().TraversePostOrder();
  auto lhs_it = lhs_order.begin();
  auto rhs_it = rhs_order.begin();
  while (lhs_it != lhs_order.end() && rhs_it != rhs_order.end()) {
    auto [lhs_node, lhs_edge] = *lhs_it;
    auto [rhs_node, rhs_edge] = *rhs_it;
    assert_equal(lhs_node.GetId().value, rhs_node.GetId().value, "Node id");
    assert_equal(lhs_edge.GetId().value, rhs_edge.GetId().value, "Edge id");
    assert_equal(lhs_node.GetCladesCount(), rhs_node.GetCladesCount(), "Clades count");
    ++lhs_it;
    ++rhs_it;
  }
}

static void test_loading_tree(std::string_view path, std::string_view refseq_path) {
  MADAG tree0 = LoadTreeFromProtobuf(path, LoadReferenceSequence(refseq_path));
  StoreTreeToProtobuf(tree0, "temp_tree.pb");
  MADAG tree1 =
      LoadTreeFromProtobuf("temp_tree.pb", LoadReferenceSequence(refseq_path));
  AssertDAGsEqual(tree0, tree1);
}

static void test_loading_dag(std::string_view path) {
  MADAG tree0 = LoadDAGFromProtobuf(path);
  StoreDAGToProtobuf(tree0, "temp_tree.pb");
  MADAG tree1 = LoadDAGFromProtobuf("temp_tree.pb");
  AssertDAGsEqual(tree0, tree1);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] {
                test_loading_tree("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
                                  "data/20D_from_fasta/refseq.fasta.gz");
              },
              "Loading: tree 20D_from_fasta"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_loading_dag("data/20D_from_fasta/full_dag.pb.gz"); },
              "Loading: dag 20D_from_fasta"});
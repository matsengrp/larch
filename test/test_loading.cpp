#include <cstddef>
#include "test_common.hpp"

#include "dag_loader.hpp"
#include "subtree_weight.hpp"
#include "parsimony_score.hpp"

#include "src/matOptimize/mutation_annotated_tree.hpp"

namespace MAT = Mutation_Annotated_Tree;
void fill_static_reference_sequence(std::string_view);
MAT::Tree mat_from_dag(const MADAG&);
MADAG build_madag_from_mat(const MAT::Tree&, std::string_view);

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
    // assert_equal(lhs_node.GetId().value, rhs_node.GetId().value, "Node id");
    // assert_equal(lhs_edge.GetId().value, rhs_edge.GetId().value, "Edge id");
    assert_equal(lhs_node.GetCladesCount(), rhs_node.GetCladesCount(), "Clades count");
    const EdgeMutations& lhs_muts = lhs.GetEdgeMutations().at(lhs_edge.GetId().value);
    const EdgeMutations& rhs_muts = rhs.GetEdgeMutations().at(rhs_edge.GetId().value);
    assert_equal(lhs_muts.size(), rhs_muts.size(), "Edge mutations count");
    for (size_t i = 0; i < lhs_muts.size(); ++i) {
      assert_equal(lhs_muts.size(), rhs_muts.size(), "Edge mutation");
    }
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

  fill_static_reference_sequence(tree0.GetReferenceSequence());
  MADAG tree2 = build_madag_from_mat(mat_from_dag(tree0), tree0.GetReferenceSequence());
  AssertDAGsEqual(tree0, tree2);

  MADAG tree3 = build_madag_from_mat(mat_from_dag(tree1), tree1.GetReferenceSequence());
  AssertDAGsEqual(tree0, tree3);
}

static void test_loading_dag(std::string_view path) {
  MADAG tree0 = LoadDAGFromProtobuf(path);
  StoreDAGToProtobuf(tree0, "temp_tree.pb");
  MADAG tree1 = LoadDAGFromProtobuf("temp_tree.pb");
  AssertDAGsEqual(tree0, tree1);

  SubtreeWeight<ParsimonyScore> weight(tree0);
  MADAG sampled0 = weight.SampleTree({}).first;

  fill_static_reference_sequence(sampled0.GetReferenceSequence());
  MADAG sampled1 =
      build_madag_from_mat(mat_from_dag(sampled0), sampled0.GetReferenceSequence());
  AssertDAGsEqual(sampled0, sampled1);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] {
                test_loading_tree("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
                                  "data/20D_from_fasta/refseq.txt.gz");
              },
              "Loading: tree 20D_from_fasta"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_loading_dag("data/20D_from_fasta/full_dag.pb.gz"); },
              "Loading: dag 20D_from_fasta"});

[[maybe_unused]] static const auto test_added2 =
    add_test({[] {
                test_loading_tree("data/startmat/startmat_no_ancestral.pb.gz",
                                  "data/startmat/refseq.txt.gz");
              },
              "Loading: tree startmat"});
#include <cstddef>
#include "test_common.hpp"

#include "larch/dag_loader.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"

#include "larch/usher_glue.hpp"

template <typename LHS, typename RHS>
static void AssertNodesEqual(LHS lhs_node, RHS rhs_node) {
  // assert_equal(lhs_node.GetId().value, rhs_node.GetId().value, "Node id");
  assert_equal(lhs_node.GetCladesCount(), rhs_node.GetCladesCount(), "Clades count");
}

template <typename LHS, typename RHS>
static void AssertEdgesEqual(LHS lhs_edge, RHS rhs_edge) {
  AssertNodesEqual(lhs_edge.GetParent(), rhs_edge.GetParent());
  AssertNodesEqual(lhs_edge.GetChild(), rhs_edge.GetChild());
  // assert_equal(lhs_edge.GetId().value, rhs_edge.GetId().value, "Edge id");

  const EdgeMutations& lhs_muts = lhs_edge.GetEdgeMutations();
  const EdgeMutations& rhs_muts = rhs_edge.GetEdgeMutations();
  assert_equal(lhs_muts.size(), rhs_muts.size(), "Edge mutations count");
  for (size_t i = 0; i < lhs_muts.size(); ++i) {
    assert_equal(lhs_muts.size(), rhs_muts.size(), "Edge mutation");
  }
}

template <typename LHS, typename RHS>
static void AssertDAGsEqual(LHS lhs, RHS rhs) {
  assert_equal(lhs.GetNodesCount(), rhs.GetNodesCount(), "Nodes count");
  assert_equal(lhs.GetEdgesCount(), rhs.GetEdgesCount(), "Edges count");

  AssertEdgesEqual(lhs.GetRoot().GetFirstChild(), rhs.GetRoot().GetFirstChild());
}

#ifdef USE_USHER
static void test_loading_tree(std::string_view path, std::string_view refseq_path) {
  MADAGStorage tree0 = LoadTreeFromProtobuf(path, LoadReferenceSequence(refseq_path));
  StoreTreeToProtobuf(tree0.View(), "temp_tree.pb");
  MADAGStorage tree1 =
      LoadTreeFromProtobuf("temp_tree.pb", LoadReferenceSequence(refseq_path));
  AssertDAGsEqual(tree0.View(), tree1.View());

  fill_static_reference_sequence(tree0.View().GetReferenceSequence());

  auto tree2 = AddMATConversion(MADAGStorage{{}});
  MAT::Tree mat0;
  AddMATConversion(tree0.View()).View().BuildMAT(mat0);
  tree2.View().BuildFromMAT(mat0, tree0.View().GetReferenceSequence());

  AssertDAGsEqual(tree0.View(), tree2.View());

  auto tree3 = AddMATConversion(MADAGStorage{{}});
  MAT::Tree mat1;
  AddMATConversion(tree1.View()).View().BuildMAT(mat1);
  tree3.View().BuildFromMAT(mat1, tree1.View().GetReferenceSequence());

  AssertDAGsEqual(tree0.View(), tree3.View());
}
#endif

static void test_loading_dag(std::string_view path) {
  MADAGStorage tree0 = LoadDAGFromProtobuf(path);
  StoreDAGToProtobuf(tree0.View(), "temp_tree.pb");
  MADAGStorage tree1 = LoadDAGFromProtobuf("temp_tree.pb");
  AssertDAGsEqual(tree0.View(), tree1.View());

  SubtreeWeight<ParsimonyScore, MADAG> weight(tree0.View());
  auto sampled0 = weight.SampleTree({});

  fill_static_reference_sequence(sampled0.View().GetReferenceSequence());
  auto sampled1 = AddMATConversion(MADAGStorage{{}});
  MAT::Tree mat0;
  AddMATConversion(sampled0.View()).View().BuildMAT(mat0);
  sampled1.View().BuildFromMAT(mat0, sampled0.View().GetReferenceSequence());
  AssertDAGsEqual(sampled0.View(), sampled1.View());
}

#ifdef USE_USHER
[[maybe_unused]] static const auto test_added0 =
    add_test({[] {
                test_loading_tree("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
                                  "data/20D_from_fasta/refseq.txt.gz");
              },
              "Loading: tree 20D_from_fasta"});
#endif

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_loading_dag("data/20D_from_fasta/full_dag.pb.gz"); },
              "Loading: dag 20D_from_fasta"});

#ifdef USE_USHER
[[maybe_unused]] static const auto test_added2 =
    add_test({[] {
                test_loading_tree("data/startmat/startmat_no_ancestral.pb.gz",
                                  "data/startmat/refseq.txt.gz");
              },
              "Loading: tree startmat"});
#endif

// [[maybe_unused]] static const auto test_added3 =
//     add_test({[] {
//                 test_loading_tree("data/20B/20B_start_tree_no_ancestral.pb.gz",
//                                   "data/20B/ref_seq_noancestral.txt.gz");
//               },
//               "Loading: tree 20B"});

// [[maybe_unused]] static const auto test_added4 =
//     add_test({[] {
//                 test_loading_tree("data/20C/20C_start_tree_no_ancestral.pb.gz",
//                                   "data/20C/ref_seq_noancestral.txt.gz");
//               },
//               "Loading: tree 20C"});

#ifdef USE_USHER
[[maybe_unused]] static const auto test_added5 =
    add_test({[] {
                test_loading_tree("data/AY.103/AY.103_start_tree_no_ancestral.pb.gz",
                                  "data/AY.103/ref_seq_noancestral.txt.gz");
              },
              "Loading: tree AY.103"});
#endif

// [[maybe_unused]] static const auto test_added6 = add_test(
//     {[] {
//        test_loading_tree("data/B.1.1.529/B.1.1.529_start_tree_no_ancestral.pb.gz",
//                          "data/B.1.1.529/ref_seq_noancestral.txt.gz");
//      },
//      "Loading: tree B.1.1.529"});

#include <cstddef>
#include "test_common.hpp"

#include "larch/dag_loader.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"

#include "larch/usher_glue.hpp"

template <typename LHS, typename RHS>
static void AssertNodesEqual(LHS lhs_node, RHS rhs_node) {
  // TestAssert(lhs_node.GetId().value == rhs_node.GetId().value);
  TestAssert(lhs_node.GetCladesCount() == rhs_node.GetCladesCount());
}

template <typename LHS, typename RHS>
static void AssertEdgesEqual(LHS lhs_edge, RHS rhs_edge) {
  AssertNodesEqual(lhs_edge.GetParent(), rhs_edge.GetParent());
  AssertNodesEqual(lhs_edge.GetChild(), rhs_edge.GetChild());
  // TestAssert(lhs_edge.GetId().value == rhs_edge.GetId().value);

  const EdgeMutations& lhs_muts = lhs_edge.GetEdgeMutations();
  const EdgeMutations& rhs_muts = rhs_edge.GetEdgeMutations();
  TestAssert(lhs_muts.size() == rhs_muts.size());
  for (size_t i = 0; i < lhs_muts.size(); ++i) {
    TestAssert(lhs_muts.size() == rhs_muts.size());
  }
}

template <typename LHS, typename RHS>
static void AssertDAGsEqual(LHS lhs, RHS rhs) {
  TestAssert(lhs.GetNodesCount() == rhs.GetNodesCount());
  TestAssert(lhs.GetEdgesCount() == rhs.GetEdgesCount());

  AssertEdgesEqual(lhs.GetRoot().GetFirstChild(), rhs.GetRoot().GetFirstChild());
}

static void test_loading_tree(std::string_view path, std::string_view refseq_path) {
  MADAGStorage tree0 = LoadTreeFromProtobuf(path, LoadReferenceSequence(refseq_path));
  StoreTreeToProtobuf(tree0.View(), test_output_folder + "/temp_tree.pb");
  MADAGStorage tree1 = LoadTreeFromProtobuf(test_output_folder + "/temp_tree.pb",
                                            LoadReferenceSequence(refseq_path));
  AssertDAGsEqual(tree0.View(), tree1.View());

  fill_static_reference_sequence(tree0.View().GetReferenceSequence());

  auto tree2 = AddMATConversion(MADAGStorage<>::EmptyDefault());
  MAT::Tree mat0;
  AddMATConversion(tree0.View()).View().BuildMAT(mat0);
  tree2.View().BuildFromMAT(mat0, tree0.View().GetReferenceSequence());

  AssertDAGsEqual(tree0.View(), tree2.View());

  auto tree3 = AddMATConversion(MADAGStorage<>::EmptyDefault());
  MAT::Tree mat1;
  AddMATConversion(tree1.View()).View().BuildMAT(mat1);
  tree3.View().BuildFromMAT(mat1, tree1.View().GetReferenceSequence());

  AssertDAGsEqual(tree0.View(), tree3.View());
}

static void test_loading_dag(std::string_view path) {
  MADAGStorage<> tree0 = LoadDAGFromProtobuf(path);
  StoreDAGToProtobuf(tree0.View(), test_output_folder + "/temp_tree.pb");
  MADAGStorage<> tree1 = LoadDAGFromProtobuf(test_output_folder + "/temp_tree.pb");
  AssertDAGsEqual(tree0.View(), tree1.View());

  SubtreeWeight<ParsimonyScore, MADAG> weight(tree0.View());
  auto sampled0 = weight.SampleTree({});

  fill_static_reference_sequence(sampled0.View().GetReferenceSequence());
  auto sampled1 = AddMATConversion(MADAGStorage<>::EmptyDefault());
  MAT::Tree mat0;
  AddMATConversion(sampled0.View()).View().BuildMAT(mat0);
  sampled1.View().BuildFromMAT(mat0, sampled0.View().GetReferenceSequence());
  AssertDAGsEqual(sampled0.View(), sampled1.View());
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] {
                test_loading_tree("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
                                  "data/20D_from_fasta/refseq.txt.gz");
              },
              "Load protobuf: tree 20D_from_fasta"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_loading_dag("data/20D_from_fasta/full_dag.pb.gz"); },
              "Load protobuf: dag 20D_from_fasta"});

[[maybe_unused]] static const auto test_added2 =
    add_test({[] {
                test_loading_tree("data/startmat/startmat_no_ancestral.pb.gz",
                                  "data/startmat/refseq.txt.gz");
              },
              "Load protobuf: tree startmat",
              {"slow"}});

// [[maybe_unused]] static const auto test_added3 =
//     add_test({[] {
//                 test_loading_tree("data/20B/20B_start_tree_no_ancestral.pb.gz",
//                                   "data/20B/ref_seq_noancestral.txt.gz");
//               },
//               "Load protobuf: tree 20B"});

// [[maybe_unused]] static const auto test_added4 =
//     add_test({[] {
//                 test_loading_tree("data/20C/20C_start_tree_no_ancestral.pb.gz",
//                                   "data/20C/ref_seq_noancestral.txt.gz");
//               },
//               "Load protobuf: tree 20C"});

[[maybe_unused]] static const auto test_added5 =
    add_test({[] {
                test_loading_tree("data/AY.103/AY.103_start_tree_no_ancestral.pb.gz",
                                  "data/AY.103/ref_seq_noancestral.txt.gz");
              },
              "Load protobuf: tree AY.103",
              {"slow"}});

// [[maybe_unused]] static const auto test_added6 = add_test(
//     {[] {
//        test_loading_tree("data/B.1.1.529/B.1.1.529_start_tree_no_ancestral.pb.gz",
//                          "data/B.1.1.529/ref_seq_noancestral.txt.gz");
//      },
//      "Load protobuf: tree B.1.1.529"});

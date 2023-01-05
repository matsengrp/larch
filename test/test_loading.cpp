#include <cstddef>
#include "test_common.hpp"

#include "larch/dag_loader.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"

#include "larch/usher_glue.hpp"

using Node = MADAG::NodeView;
using Edge = MADAG::EdgeView;

static void AssertEqual(Node lhs_node, Node rhs_node) {
  // assert_equal(lhs_node.GetId().value, rhs_node.GetId().value, "Node id");
  assert_equal(lhs_node.GetCladesCount(), rhs_node.GetCladesCount(), "Clades count");
}

static void AssertEqual(Edge lhs_edge, Edge rhs_edge) {
  AssertEqual(lhs_edge.GetParent(), rhs_edge.GetParent());
  AssertEqual(lhs_edge.GetChild(), rhs_edge.GetChild());
  // assert_equal(lhs_edge.GetId().value, rhs_edge.GetId().value, "Edge id");

  const EdgeMutations& lhs_muts = lhs_edge.GetEdgeMutations();
  const EdgeMutations& rhs_muts = rhs_edge.GetEdgeMutations();
  assert_equal(lhs_muts.size(), rhs_muts.size(), "Edge mutations count");
  for (size_t i = 0; i < lhs_muts.size(); ++i) {
    assert_equal(lhs_muts.size(), rhs_muts.size(), "Edge mutation");
  }
}

static void AssertDAGsEqual(MADAG lhs, MADAG rhs) {
  assert_equal(lhs.GetNodesCount(), rhs.GetNodesCount(), "Nodes count");
  assert_equal(lhs.GetEdgesCount(), rhs.GetEdgesCount(), "Edges count");

  AssertEqual(lhs.GetRoot().GetFirstChild(), rhs.GetRoot().GetFirstChild());
}

static void test_loading_tree(std::string_view path, std::string_view refseq_path) {
  MADAGStorage tree0 = LoadTreeFromProtobuf(path, LoadReferenceSequence(refseq_path));
  StoreTreeToProtobuf(tree0.View(), "temp_tree.pb");
  MADAGStorage tree1 =
      LoadTreeFromProtobuf("temp_tree.pb", LoadReferenceSequence(refseq_path));
  AssertDAGsEqual(tree0.View(), tree1.View());

  fill_static_reference_sequence(tree0.View().GetReferenceSequence());
  MADAGStorage tree2 = build_madag_from_mat(mat_from_dag(tree0.View()),
                                            tree0.View().GetReferenceSequence());
  AssertDAGsEqual(tree0.View(), tree2.View());

  MADAGStorage tree3 = build_madag_from_mat(mat_from_dag(tree1.View()),
                                            tree1.View().GetReferenceSequence());
  AssertDAGsEqual(tree0.View(), tree3.View());
}

static void test_loading_dag(std::string_view path) {
  MADAGStorage tree0 = LoadDAGFromProtobuf(path);
  StoreDAGToProtobuf(tree0.View(), "temp_tree.pb");
  MADAGStorage tree1 = LoadDAGFromProtobuf("temp_tree.pb");
  AssertDAGsEqual(tree0.View(), tree1.View());

  SubtreeWeight<ParsimonyScore, MADAG> weight(tree0.View());
  MADAGStorage sampled0 = weight.SampleTree({}).first;

  fill_static_reference_sequence(sampled0.View().GetReferenceSequence());
  MADAGStorage sampled1 = build_madag_from_mat(mat_from_dag(sampled0.View()),
                                               sampled0.View().GetReferenceSequence());
  AssertDAGsEqual(sampled0.View(), sampled1.View());
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

[[maybe_unused]] static const auto test_added3 =
    add_test({[] {
                test_loading_tree("data/20B/20B_start_tree_no_ancestral.pb.gz",
                                  "data/20B/ref_seq_noancestral.txt.gz");
              },
              "Loading: tree 20B"});

[[maybe_unused]] static const auto test_added4 =
    add_test({[] {
                test_loading_tree("data/20C/20C_start_tree_no_ancestral.pb.gz",
                                  "data/20C/ref_seq_noancestral.txt.gz");
              },
              "Loading: tree 20C"});

[[maybe_unused]] static const auto test_added5 =
    add_test({[] {
                test_loading_tree("data/AY.103/AY.103_start_tree_no_ancestral.pb.gz",
                                  "data/AY.103/ref_seq_noancestral.txt.gz");
              },
              "Loading: tree AY.103"});

[[maybe_unused]] static const auto test_added6 = add_test(
    {[] {
       test_loading_tree("data/B.1.1.529/B.1.1.529_start_tree_no_ancestral.pb.gz",
                         "data/B.1.1.529/ref_seq_noancestral.txt.gz");
     },
     "Loading: tree B.1.1.529"});

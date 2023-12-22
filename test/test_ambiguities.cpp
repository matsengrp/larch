#include "test_common.hpp"
#include "sample_dag.hpp"
#include "larch/dag_loader.hpp"

[[maybe_unused]] static auto BuildNodeSequenceMap(MADAGStorage<> &dag_storage,
                                                  bool include_nonleaf_nodes = false) {
  NodeSeqMap node_seq_map;
  auto dag = dag_storage.View();
  auto ref_seq = dag.GetReferenceSequence();
  for (auto node : dag.GetNodes()) {
    if (include_nonleaf_nodes || node.IsLeaf()) {
      node_seq_map[node.GetId()] = node.GetCompactGenome().ToSequence(ref_seq);
    }
  }
  return node_seq_map;
}

[[maybe_unused]] static auto DAGMutationsAreValid(MADAGStorage<> &dag_storage) {
  using Edge = MutableMADAG::EdgeView;
  auto dag = dag_storage.View();

  auto EdgeIsValid = [&dag](const Edge edge) {
    const auto ref_seq = dag.GetReferenceSequence();
    auto parent = dag.Get(edge.GetParent());
    auto child = dag.Get(edge.GetChild());
    const auto &edge_muts = edge.GetEdgeMutations();
    const auto &parent_cg = parent.GetCompactGenome();
    const auto &child_cg = child.GetCompactGenome();
    auto test_muts = CompactGenome::ToEdgeMutations(ref_seq, parent_cg, child_cg);
    return (edge_muts == test_muts);
  };

  for (auto edge : dag.GetEdges()) {
    if (!EdgeIsValid(edge)) return false;
  }
  return true;
}

[[maybe_unused]] static auto VerifyCompactGenomesCompatibleWithLeaves(
    MADAGStorage<> &dag_storage, NodeSeqMap &truth_leaf_seq_map) {
  auto dag = dag_storage.View();
  auto dag_leaf_seq_map = BuildNodeSequenceMap(dag_storage, false);
  for (auto node : dag.GetLeafs()) {
    if (dag_leaf_seq_map[node.GetId()] != truth_leaf_seq_map[node.GetId()]) {
      return false;
    }
  }
  return true;
}

[[maybe_unused]] static void WriteDAGToFile(MADAGStorage<> &dag_storage,
                                            const std::string &output_filename) {
  std::ofstream os;
  auto dag = dag_storage.View();
  os.open(output_filename);
  MADAGToDOT(dag, os);
  os.close();
}

[[maybe_unused]] void test_compare_ambiguities() {
  auto amb_dag_storage = make_ambiguous_sample_dag();
  auto amb_seq_map = make_sample_ambiguous_sequence_map();
  auto amb_dag = amb_dag_storage.View();

  auto unamb_dag_storage = make_unambiguous_sample_dag();
  auto unamb_seq_map = make_sample_unambiguous_sequence_map();
  auto unamb_dag = unamb_dag_storage.View();

  bool write_files = false;
  if (write_files) {
    WriteDAGToFile(amb_dag_storage, "_ignore/amb_dag.dot");
    WriteDAGToFile(unamb_dag_storage, "_ignore/unamb_dag.dot");
    StoreDAGToProtobuf(amb_dag, "_ignore/amb_dag.pb");
    StoreDAGToProtobuf(unamb_dag, "_ignore/unamb_dag.pb");
  }

  // (0) Test checks that all edges mutations are compatible with adjacent compact
  // genomes.
  amb_dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{1}] = {'T', 'G'};
  assert_false(DAGMutationsAreValid(amb_dag_storage),
               "Test_0a: DAG EdgeMutations incorrectly found to be compatible with "
               "Compact Genomes.");
  amb_dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{1}] = {'G', 'A'};
  assert_false(DAGMutationsAreValid(amb_dag_storage),
               "Test_0b: DAG EdgeMutations incorrectly found to be compatible with "
               "Compact Genomes.");
  amb_dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{1}] = {'T', 'A'};
  assert_true(DAGMutationsAreValid(amb_dag_storage),
              "Test_0c: DAG EdgeMutations are not compatible with Compact Genomes.");
  assert_true(DAGMutationsAreValid(unamb_dag_storage),
              "Test_0d: DAG EdgeMutations are not compatible with Compact Genomes.");

  // (1) Test checks for ambiguity.
  for (auto base : {'A', 'C', 'G', 'T'}) {
    assert_false(MutationBase{base}.IsAmbiguous(),
                 "Test_0a: MutationBase found incorrectly ambiguous.");
  }
  assert_true(MutationBase{'N'}.IsAmbiguous(),
              "Test_1b: MutationBase found incorrectly unambiguous.");
  assert_true(MutationBase{{0, 1, 0, 1}}.IsAmbiguous(),
              "Test_1c: MutationBase found incorrectly unambiguous.");

  // (2) Test ambiguous comparisons of bases.
  for (auto base : {'A', 'C', 'G', 'T'}) {
    assert_true(MutationBase{base}.IsCompatible(MutationBase{'N'}),
                "Test_2a: AmbiguousCompare found incorrectly incompatible.");
    assert_true(MutationBase{'N'}.IsCompatible(MutationBase{base}),
                "Test_2b: AmbiguousCompare found incorrectly incompatible.");
  }
  assert_true(MutationBase{{1, 1, 1, 1}}.IsCompatible(MutationBase{{1, 1, 1, 1}}),
              "Test_2c: AmbiguousCompare found incorrectly incompatible.");
  assert_false(MutationBase{{0, 0, 1, 1}}.IsCompatible(MutationBase{{1, 1, 0, 0}}),
               "Test_2d: AmbiguousCompare found incorrectly compatible.");

  assert_true(MutationBase{{0, 1, 1, 0}}.GetFirstBase() == MutationBase{{0, 1, 0, 0}},
              "Test_2e: GetFirstBase does not return correct value.");
  assert_true(MutationBase{{0, 1, 1, 0}}.GetFirstCommonBase(
                  MutationBase{{0, 0, 1, 1}}) == MutationBase{{0, 0, 1, 0}},
              "Test_2f: GetFirstCommonBase does not return correct value.");
  assert_true(MutationBase{{1, 1, 1, 0}}.GetCommonBases(MutationBase{{0, 1, 1, 1}}) ==
                  MutationBase{{0, 1, 1, 0}},
              "Test_2g: GetCommonBases does not return correct value.");

  // (3) Test that ambiguous leaves are compatible with unambiguous leaves.
  for (auto node : amb_dag.GetLeafs()) {
    auto &cg_1 = amb_dag.Get(node.GetId()).GetCompactGenome();
    auto &cg_2 = unamb_dag.Get(node.GetId()).GetCompactGenome();
    bool cgs_compatible = cg_1.IsCompatible(cg_2, unamb_dag.GetReferenceSequence());
    cgs_compatible &= cg_2.IsCompatible(cg_1, unamb_dag.GetReferenceSequence());
    assert_true(cgs_compatible,
                "Test_3a: Ambiguous Compact Genomes incorrectly found incompatible.");
  }
  auto &cg_1 = unamb_dag.Get(NodeId{1}).GetCompactGenome();
  auto &cg_2 = unamb_dag.Get(NodeId{2}).GetCompactGenome();
  bool cgs_compatible = cg_1.IsCompatible(cg_2, unamb_dag.GetReferenceSequence());
  cgs_compatible &= cg_1.IsCompatible(cg_2, unamb_dag.GetReferenceSequence());
  assert_false(cgs_compatible,
               "Test_3b: Ambiguous Compact Genomes incorrectly found compatible.");

  // (4) Test that edge ambiguous leaves don't form edge mutations.
  std::vector<EdgeId> amb_edge_ids{{{2}, {4}}};
  for (auto edge_id : amb_edge_ids) {
    auto &edge_mut_1 = unamb_dag.Get(edge_id).GetEdgeMutations();
    auto &edge_mut_2 = amb_dag.Get(edge_id).GetEdgeMutations();
    assert_true(edge_mut_1.size() > edge_mut_2.size(),
                "Test_4a: Number of Ambiguous Edge Mutations incorrectly >= Unambigous "
                "Edge Mutations.");
  }
  std::vector<EdgeId> unamb_edge_ids{{{1}, {3}}};
  for (auto edge_id : unamb_edge_ids) {
    auto &edge_mut_1 = unamb_dag.Get(edge_id).GetEdgeMutations();
    auto &edge_mut_2 = amb_dag.Get(edge_id).GetEdgeMutations();
    assert_true(edge_mut_1.size() == edge_mut_2.size(),
                "Test_4b: Number of Ambiguous Edge Mutations with no ambiguities "
                "incorrectly != Unambigous Edge Mutations.");
  }

  // (5) Test known DAG edge mutation counts.
  // 'TGG' -> 'TAG'/'TNN'
  assert_true(
      unamb_dag.Get(EdgeId{2}).GetEdgeMutations().size() == 1,
      "Test_5a: Number of Unambiguous Edge Mutations does not match known number.");
  assert_true(
      amb_dag.Get(EdgeId{2}).GetEdgeMutations().size() == 0,
      "Test_5b: Number of Ambiguous Edge Mutations does not match known number.");
  // 'CTC' -> 'ACG'/'ANG'
  assert_true(
      unamb_dag.Get(EdgeId{4}).GetEdgeMutations().size() == 3,
      "Test_5c: Number of Unambiguous Edge Mutations does not match known number.");
  assert_true(
      amb_dag.Get(EdgeId{4}).GetEdgeMutations().size() == 2,
      "Test_5d: Number of Ambiguous Edge Mutations does not match known number.");

  // (6) Test that recomputing edge mutations and compact genomes does not alter
  // leaves.
  amb_dag.RecomputeCompactGenomes(false);
  assert_true(VerifyCompactGenomesCompatibleWithLeaves(amb_dag_storage, amb_seq_map),
              "Test_6a: RecomputeCompactGenomes incorrectly altered leaf CGs.");
  unamb_dag.RecomputeCompactGenomes(false);
  assert_true(
      VerifyCompactGenomesCompatibleWithLeaves(unamb_dag_storage, unamb_seq_map),
      "Test_6b: RecomputeCompactGenomes incorrectly altered leaf CGs.");
  amb_dag.RecomputeCompactGenomes(true);
  assert_false(
      VerifyCompactGenomesCompatibleWithLeaves(amb_dag_storage, amb_seq_map),
      "Test_6c: RecomputeCompactGenomes incorrectly left leaf CGs unaltered .");
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[]() { test_compare_ambiguities(); }, "Ambiguities Test"});

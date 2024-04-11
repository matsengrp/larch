#include "test_common.hpp"
#include "sample_dag.hpp"
#include "larch/dag_loader.hpp"

[[maybe_unused]] static auto build_node_sequence_map(
    MADAGStorage<> &dag_storage, bool include_nonleaf_nodes = false) {
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

[[maybe_unused]] static auto verify_dag_mutations_are_valid(
    MADAGStorage<> &dag_storage) {
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

[[maybe_unused]] static auto verify_compact_genomes_compatible_with_leaves(
    MADAGStorage<> &dag_storage, NodeSeqMap &truth_leaf_seq_map,
    bool do_print_failure = true) {
  auto dag = dag_storage.View();
  auto dag_leaf_seq_map = build_node_sequence_map(dag_storage, false);
  for (auto node : dag.GetLeafs()) {
    if (dag_leaf_seq_map[node.GetId()] != truth_leaf_seq_map[node.GetId()]) {
      if (do_print_failure) {
        std::cout << "Failed_at [" << node.GetId()
                  << "]: TEST:" << truth_leaf_seq_map[node.GetId()]
                  << " vs TRUTH:" << dag_leaf_seq_map[node.GetId()] << std::endl;
      }
      return false;
    }
  }
  return true;
}

[[maybe_unused]] static void write_dag_to_file(MADAGStorage<> &dag_storage,
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
    write_dag_to_file(amb_dag_storage, "_ignore/amb_dag.dot");
    write_dag_to_file(unamb_dag_storage, "_ignore/unamb_dag.dot");
    StoreDAGToProtobuf(amb_dag, "_ignore/amb_dag.pb");
    StoreDAGToProtobuf(unamb_dag, "_ignore/unamb_dag.pb");
  }

  // (0) Test spot checks that edges mutations are compatible with adjacent compact
  // genomes.
  // Initial edge mutation: T->A.
  // Set to T->G.
  amb_dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{1}] = {'T', 'G'};
  TestAssert(not(verify_dag_mutations_are_valid(amb_dag_storage)) &&
             "Test_0a: DAG EdgeMutations incorrectly found to be compatible with "
             "Compact Genomes.");
  // Set to G->A.
  amb_dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{1}] = {'G', 'A'};
  TestAssert(not(verify_dag_mutations_are_valid(amb_dag_storage)) &&
             "Test_0b: DAG EdgeMutations incorrectly found to be compatible with "
             "Compact Genomes.");
  // Reset back to T->A.
  amb_dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{1}] = {'T', 'A'};
  TestAssert((verify_dag_mutations_are_valid(amb_dag_storage)) &&
             "Test_0c: DAG EdgeMutations are not compatible with ambiguous Compact "
             "Genomes.");
  TestAssert((verify_dag_mutations_are_valid(unamb_dag_storage)) &&
             "Test_0d: DAG EdgeMutations are not compatible with unambiguous Compact "
             "Genomes.");

  // (1) Test checks for ambiguity.
  for (auto base : {'A', 'C', 'G', 'T'}) {
    TestAssert(not(MutationBase{base}.IsAmbiguous()) &&
               "Test_0a: MutationBase found incorrectly ambiguous.");
  }
  TestAssert((MutationBase{'N'}.IsAmbiguous()) &&
             "Test_1b: MutationBase found incorrectly unambiguous.");
  TestAssert((MutationBase{{0, 1, 0, 1}}.IsAmbiguous()) &&
             "Test_1c: MutationBase found incorrectly unambiguous.");

  // (2) Test ambiguous comparisons of bases.
  for (auto base : {'A', 'C', 'G', 'T'}) {
    TestAssert((MutationBase{base}.IsCompatible(MutationBase{'N'})) &&
               "Test_2a: AmbiguousCompare found incorrectly incompatible.");
    TestAssert((MutationBase{'N'}.IsCompatible(MutationBase{base})) &&
               "Test_2b: AmbiguousCompare found incorrectly incompatible.");
  }
  TestAssert((MutationBase{{1, 1, 1, 1}}.IsCompatible(MutationBase{{1, 1, 1, 1}})) &&
             "Test_2c: AmbiguousCompare found incorrectly incompatible.");
  TestAssert(not(MutationBase{{0, 0, 1, 1}}.IsCompatible(MutationBase{{1, 1, 0, 0}})) &&
             "Test_2d: AmbiguousCompare found incorrectly compatible.");

  TestAssert(
      (MutationBase{{0, 1, 1, 0}}.GetFirstBase() == MutationBase{{0, 1, 0, 0}}) &&
      "Test_2e: GetFirstBase does not return correct value.");
  TestAssert((MutationBase{{0, 1, 1, 0}}.GetFirstCommonBase(
                  MutationBase{{0, 0, 1, 1}}) == MutationBase{{0, 0, 1, 0}}) &&
             "Test_2f: GetFirstCommonBase does not return correct value.");
  TestAssert((MutationBase{{1, 1, 1, 0}}.GetCommonBases(MutationBase{{0, 1, 1, 1}}) ==
              MutationBase{{0, 1, 1, 0}}) &&
             "Test_2g: GetCommonBases does not return correct value.");

  // (3) Test that ambiguous leaves are compatible with unambiguous leaves.

  for (auto node : amb_dag.GetLeafs()) {
    auto &cg_1 = amb_dag.Get(node.GetId()).GetCompactGenome();
    auto &cg_2 = unamb_dag.Get(node.GetId()).GetCompactGenome();
    bool cgs_compatible = cg_1.IsCompatible(cg_2, unamb_dag.GetReferenceSequence());
    cgs_compatible &= cg_2.IsCompatible(cg_1, unamb_dag.GetReferenceSequence());
    TestAssert((cgs_compatible) &&
               "Test_3a: Ambiguous Compact Genomes incorrectly found incompatible.");
  }
  auto &cg_1 = unamb_dag.Get(NodeId{1}).GetCompactGenome();
  auto &cg_2 = unamb_dag.Get(NodeId{2}).GetCompactGenome();
  bool cgs_compatible = cg_1.IsCompatible(cg_2, unamb_dag.GetReferenceSequence());
  cgs_compatible &= cg_1.IsCompatible(cg_2, unamb_dag.GetReferenceSequence());
  TestAssert(not(cgs_compatible) &&
             "Test_3b: Ambiguous Compact Genomes incorrectly found compatible.");

  // (4) Test that edge ambiguous leaves don't form edge mutations.
  std::vector<EdgeId> amb_edge_ids{{{2}, {4}}};
  for (auto edge_id : amb_edge_ids) {
    auto &edge_mut_1 = unamb_dag.Get(edge_id).GetEdgeMutations();
    auto &edge_mut_2 = amb_dag.Get(edge_id).GetEdgeMutations();
    TestAssert((edge_mut_1.size() > edge_mut_2.size()) &&
               "Test_4a: Number of Ambiguous Edge Mutations incorrectly >= Unambigous "
               "Edge Mutations.");
  }
  std::vector<EdgeId> unamb_edge_ids{{{1}, {3}}};
  for (auto edge_id : unamb_edge_ids) {
    auto &edge_mut_1 = unamb_dag.Get(edge_id).GetEdgeMutations();
    auto &edge_mut_2 = amb_dag.Get(edge_id).GetEdgeMutations();
    TestAssert((edge_mut_1.size() == edge_mut_2.size()) &&
               "Test_4b: Number of Ambiguous Edge Mutations with no ambiguities "
               "incorrectly != Unambigous Edge Mutations.");
  }

  // (5) Test known DAG edge mutation counts.
  // 'TGG' -> 'TAG'/'TNN'
  TestAssert(
      (unamb_dag.Get(EdgeId{2}).GetEdgeMutations().size() == 1) &&
      "Test_5a: Number of Unambiguous Edge Mutations does not match known number.");
  TestAssert(
      (amb_dag.Get(EdgeId{2}).GetEdgeMutations().size() == 0) &&
      "Test_5b: Number of Ambiguous Edge Mutations does not match known number.");
  // 'CTC' -> 'ACG'/'ANG'
  TestAssert(
      (unamb_dag.Get(EdgeId{4}).GetEdgeMutations().size() == 3) &&
      "Test_5c: Number of Unambiguous Edge Mutations does not match known number.");
  TestAssert(
      (amb_dag.Get(EdgeId{4}).GetEdgeMutations().size() == 2) &&
      "Test_5d: Number of Ambiguous Edge Mutations does not match known number.");

  // (6) Test that recomputing edge mutations and compact genomes does not alter
  // leaves.
  TestAssert(
      (verify_compact_genomes_compatible_with_leaves(amb_dag_storage, amb_seq_map)) &&
      "Test_6a: Leaf CGs altered before running RecomputeCompactGenomes.");
  amb_dag.RecomputeCompactGenomes(false);
  verify_compact_genomes_compatible_with_leaves(amb_dag_storage, amb_seq_map);
  TestAssert(
      (verify_compact_genomes_compatible_with_leaves(amb_dag_storage, amb_seq_map)) &&
      "Test_6b: RecomputeCompactGenomes incorrectly altered leaf CGs.");
  unamb_dag.RecomputeCompactGenomes(false);
  TestAssert((verify_compact_genomes_compatible_with_leaves(unamb_dag_storage,
                                                            unamb_seq_map)) &&
             "Test_6c: RecomputeCompactGenomes incorrectly altered leaf CGs.");
  amb_dag.RecomputeCompactGenomes(true);
  TestAssert(not(verify_compact_genomes_compatible_with_leaves(amb_dag_storage,
                                                               amb_seq_map, false)) &&
             "Test_6d: RecomputeCompactGenomes incorrectly left leaf CGs unaltered .");
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[]() { test_compare_ambiguities(); }, "Ambiguities Test"});

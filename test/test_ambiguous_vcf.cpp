#include "test_common.hpp"
#include "sample_dag.hpp"

#include "larch/dag_loader.hpp"
#include "larch/spr/spr_view.hpp"
#include "larch/mat_conversion.hpp"

#include "larch/usher_glue.hpp"
// #include "larch/import_vcf.hpp"

std::ostream &operator<<(std::ostream &os, const MAT::Mutation mut) {
  os << mut.get_string();
  return os;
}

[[maybe_unused]] static std::string madag_info(const MADAGStorage &dag_storage) {
  std::stringstream os;
  auto dag = dag_storage.View();
  auto ref_seq = dag.GetReferenceSequence();
  os << "ref_seq: " << ref_seq << std::endl;
  os << "=== NODES: " << dag.GetNodesCount() << " ===" << std::endl;
  os << "=== LEAFS: " << dag.GetLeafsCount() << " ===" << std::endl;
  for (const auto node : dag.GetNodes()) {
    os << "node: " << node.GetId() << " " << node.GetSampleId().value_or("<no_name>")
       << std::endl;
    os << "   children: " << node.GetChildren() << std::endl;
    os << "   seq: " << node.GetCompactGenome().ToSequence(ref_seq) << std::endl;
    os << "   muts: " << node.GetCompactGenome().ToString() << std::endl;
  }
  return os.str();
}

[[maybe_unused]] static std::string mat_info(MAT::Tree &mat) {
  std::stringstream os;
  os << "newick: " << mat.get_newick_string(true, false, false, true) << std::endl;
  os << "=== NODES ===" << std::endl;
  for (const auto node : mat.breadth_first_expansion()) {
    auto node_name = (mat.get_node_name(node->node_id) != "")
                         ? mat.get_node_name(node->node_id)
                         : "<no_name>";
    os << "node: " << node->node_id << " " << node_name << std::endl;
    os << "   children: [ ";
    for (auto child : node->children) {
      os << child->node_id << " ";
    }
    os << "]" << std::endl;
    os << "   muts: [ ";
    for (const auto mut : node->mutations) {
      os << mut.get_string() << " ";
    }
    os << "]" << std::endl;

    const auto node_from_mat = mat.get_node(node->node_id);
    if (node_from_mat == nullptr) {
      std::cout << "node " << node->node_id << " is NULL" << std::endl;
    }
  }
  return os.str();
}

[[maybe_unused]] static std::string compact_genome_data_info(
    const std::unordered_map<std::string, CompactGenomeData> &mut_map) {
  std::stringstream os;
  os << "[ ";
  for (auto &[name, cg_data] : mut_map) {
    os << "(" << name << ": ";
    for (auto &[pos, mut] : cg_data) {
      os << pos << mut << " ";
    }
    os << ") ";
  }
  os << "]";
  return os.str();
}

[[maybe_unused]] static std::string original_state_info(
    const Original_State_t &og_state) {
  std::stringstream os;
  os << "[ ";
  for (const auto &[id, mut_set] : og_state) {
    os << "( " << id << ", [ ";
    for (const auto &mut : mut_set) {
      os << mut.get_string() << " ";
    }
    os << "] ) ";
  }
  os << "]";
  return os.str();
}

[[maybe_unused]] static int original_state_compare(const Original_State_t &lhs,
                                                   const Original_State_t &rhs) {
  std::ignore = lhs;
  std::ignore = rhs;
  // for (const auto &[lhs_id, lhs_mut_set] : lhs) {
  //   const auto &[rhs_id, rhs_mut_set] = rhs[lhs_id];
  //   for (const auto &mut : mut_set) {
  //   }
  // }
  return 0;
}

bool operator==(const Original_State_t &lhs, const Original_State_t &rhs) {
  return original_state_compare(lhs, rhs) == 0;
}

[[maybe_unused]] static auto convert_madag_to_mat(const MADAGStorage &dag_storage) {
  MAT::Tree mat;
  auto dag = dag_storage.View();
  Assert(dag.IsTree());
  auto mat_conv = AddMATConversion(dag);
  mat_conv.View().BuildMAT(mat);
  return mat;
}

[[maybe_unused]] static auto load_mat_from_protobuf(std::string_view pb_path) {
  auto dag_storage = LoadDAGFromProtobuf(pb_path);
  auto mat = convert_madag_to_mat(dag_storage);
  return mat;
}

[[maybe_unused]] static void madag_label_nodes(MADAGStorage &dag_storage,
                                               bool leaves_only = true) {
  auto dag = dag_storage.View();
  for (auto node : dag.GetNodes()) {
    if (!leaves_only || node.IsLeaf()) {
      auto prefix = (node.IsLeaf() ? "x" : "y");
      prefix = (node.IsUA() ? "r" : prefix);
      auto full_str_id = prefix + std::to_string(node.GetId().value);
      if constexpr (decltype(node)::template contains_feature<Deduplicate<SampleId>>) {
        auto id_iter = dag.template AsFeature<Deduplicate<SampleId>>().AddDeduplicated(
            SampleId{full_str_id});
        node = id_iter.first;
      } else {
        node.SetSampleId({full_str_id});
      }
    }
  }
}

[[maybe_unused]] static std::string madag_to_fasta(const MADAGStorage &dag_storage,
                                                   bool use_ids, bool leaves_only,
                                                   bool include_reference) {
  std::stringstream os;
  auto dag = dag_storage.View();
  const auto ref_seq = dag.GetReferenceSequence();
  if (include_reference) {
    os << ">ref" << std::endl;
    os << ref_seq << std::endl;
  }
  for (const auto node : dag.GetNodes()) {
    if (!leaves_only || node.IsLeaf()) {
      auto name = node.GetSampleId().value_or(std::to_string(node.GetId().value));
      name = (use_ids ? std::to_string(node.GetId().value) : name);
      os << ">" << name << std::endl;
      os << node.GetCompactGenome().ToSequence(ref_seq) << std::endl;
    }
  }
  return os.str();
}

[[maybe_unused]] static void madag_to_fasta(const MADAGStorage &dag_storage,
                                            const std::string &out_path, bool use_ids,
                                            bool leaves_only, bool include_reference) {
  std::ofstream file_out;
  file_out.open(out_path);
  file_out << madag_to_fasta(dag_storage, use_ids, leaves_only, include_reference)
           << std::endl;
  file_out.close();
}

[[maybe_unused]] static void mat_apply_compact_genome_data(
    MAT::Tree &tree,
    const std::unordered_map<std::string, CompactGenomeData> &mut_map) {
  for (const auto &[name, muts] : mut_map) {
    auto node = tree.get_node(tree.node_name_to_node_idx(name));
    if (node) {
      node->mutations.clear();
      node->mutations.reserve(muts.size());
      for (const auto &[pos, base] : muts) {
        Assert(pos.value != NoId);
        MAT::Mutation mat_mut(
            "ref", static_cast<int>(pos.value), EncodeBaseMAT(base.ToChar()),
            EncodeBaseMAT(base.ToChar()), EncodeBaseMAT(base.ToChar()));
        node->mutations.push_back(mat_mut);
      }
    } else {
      std::cerr << "ERROR: Could not find sample `" << name << "` in tree."
                << std::endl;
    }
  }
}

[[maybe_unused]] static void mat_apply_vcf(MAT::Tree &tree,
                                           const std::string &vcf_path) {
  std::ignore = tree;
  std::ignore = vcf_path;
  // VCF_input(vcf_path.c_str(), tree);
}

// MAIN TESTS

[[maybe_unused]] static auto test_ambiguous_vcf() {
  // Reference Sequence (for when reading from pb)
  auto ref_seq = "GAA";
  // VCFs
  std::string amb_vcf_path = "data/test_ambiguous_vcf/amb.vcf";
  std::string unamb_vcf_path = "data/test_ambiguous_vcf/unamb.vcf";

  // Original State objects are used for testing purposes.
  Original_State_t amb_mat_truth_og, unamb_mat_truth_og;
  Original_State_t amb_mat_og, unamb_mat_og;
  Original_State_t amb_mat_cg_og, unamb_mat_cg_og;

  // Truth MAT created from ambiguous sample DAG
  auto amb_dag_storage = MakeAmbiguousSampleDAG();
  madag_label_nodes(amb_dag_storage, false);
  madag_to_fasta(amb_dag_storage, "data/test_ambiguous_vcf/amb.fasta", false, true,
                 true);
  auto amb_mat_truth = convert_madag_to_mat(amb_dag_storage);
  check_samples(amb_mat_truth.root, amb_mat_truth_og, &amb_mat_truth, true);
  // Truth MAT created from unambiguous sample DAG
  auto unamb_dag_storage = MakeUnambiguousSampleDAG();
  madag_label_nodes(unamb_dag_storage, false);
  madag_to_fasta(unamb_dag_storage, "data/test_ambiguous_vcf/unamb.fasta", false, true,
                 true);
  auto unamb_mat_truth = convert_madag_to_mat(unamb_dag_storage);
  check_samples(unamb_mat_truth.root, unamb_mat_truth_og, &unamb_mat_truth, true);

  // Create topology (empty starting point)
  auto topo_dag_storage = MakeUnambiguousSampleDAG();
  madag_label_nodes(topo_dag_storage, false);
  topo_dag_storage.View().RecomputeCompactGenomes();
  topo_dag_storage.View().RecomputeEdgeMutations();
  // Ambiguous test MATs created from empty topology.
  auto amb_mat_cg = convert_madag_to_mat(topo_dag_storage);
  // check_samples(amb_mat_cg.root, amb_mat_cg_og, &amb_mat_cg, true);
  // assert_true(amb_mat_cg_og == amb_mat_truth_og,
  //             "Ambiguous MATs incorrectly match before applying VCF data.");
  // Unambiguous test MATs created from empty topology.
  auto unamb_mat_cg = convert_madag_to_mat(topo_dag_storage);
  // check_samples(unamb_mat_cg.root, unamb_mat_og, &unamb_mat_cg, true);
  // assert_true(unamb_mat_og == unamb_mat_truth_og,
  //             "Unambiguous MATs incorrectly match before applying VCF data.");

  // Update MAT from ambiguous VCF via Compact Genome data.
  auto amb_cg = LoadCompactGenomeDataFromVCF(amb_vcf_path, ref_seq);
  std::cout << "amb_cg: " << compact_genome_data_info(amb_cg) << std::endl;
  mat_apply_compact_genome_data(amb_mat_cg, amb_cg);
  check_samples(amb_mat_cg.root, amb_mat_cg_og, &amb_mat_cg, true);
  // Update MAT from unambiguous VCF via Compact Genome data.
  auto unamb_cg = LoadCompactGenomeDataFromVCF(unamb_vcf_path, ref_seq);
  std::cout << "unamb_cg: " << compact_genome_data_info(unamb_cg) << std::endl;
  mat_apply_compact_genome_data(unamb_mat_cg, unamb_cg);
  check_samples(unamb_mat_cg.root, unamb_mat_cg_og, &unamb_mat_cg, true);

  std::cout << "amb_mat_truth: " << original_state_info(amb_mat_truth_og) << std::endl;
  std::cout << "unamb_mat_truth: " << original_state_info(unamb_mat_truth_og)
            << std::endl;

  std::cout << "amb_mat: " << original_state_info(amb_mat_og) << std::endl;
  std::cout << "unamb_mat: " << original_state_info(unamb_mat_og) << std::endl;
  std::cout << "amb_mat_cg: " << original_state_info(amb_mat_cg_og) << std::endl;
  std::cout << "unamb_mat_cg: " << original_state_info(unamb_mat_cg_og) << std::endl;

  assert_true(amb_mat_cg_og == amb_mat_truth_og,
              "Ambiguous MATs do not match after applying VCF data.");
  assert_true(unamb_mat_cg_og == unamb_mat_truth_og,
              "Unambiguous MATs do not match after applying VCF data.");
}

void test_vcf_compatible() {
  auto dag_storage = MakeSampleDAG();
  madag_label_nodes(dag_storage, false);
  MADAG dag = dag_storage.View();
  auto ref_seq = dag.GetReferenceSequence();
  auto id_to_cg_map = LoadCompactGenomeDataFromVCF(
      "data/test_ambiguous_vcf/SampleDAG_unique_ambiguous_leafs.vcf", ref_seq);

  // make sure that each of the leaves in dag match exactly one key of id_to_cg_map
  for (auto leaf_id : dag.GetLeafs()) {
    if (leaf_id.GetSampleId().has_value()) {
      Assert(id_to_cg_map.find(leaf_id.GetSampleId().value()) != id_to_cg_map.end());
    }
  }

  // make sure that each parent edge is "compatible" with the leaf nodes
  MADAGApplyCompactGenomeData(dag_storage, id_to_cg_map);
  dag_storage.View().RecomputeEdgeMutations();
  for (auto leaf : dag.GetLeafs()) {
    for (auto parent_edge : leaf.GetParents()) {
      auto leaf_cg = leaf.GetCompactGenome().Copy();
      auto parent_cg = parent_edge.GetParent().GetCompactGenome().Copy();
      auto edge_mutations_from_dag = parent_edge.GetEdgeMutations().Copy();
      ContiguousMap<MutationPosition, MutationBase> current_edge_mutations;
      for (auto posval : parent_edge.GetEdgeMutations()) {
        current_edge_mutations[posval.first] = posval.second.second;
      }
      parent_cg.ApplyChanges(current_edge_mutations);
      Assert(leaf_cg.IsCompatible(parent_cg, ref_seq));
    }
  }
}

void test_vcf_reading() {
  auto dag_storage = MakeSampleDAG();
  madag_label_nodes(dag_storage, false);
  MADAG dag = dag_storage.View();
  auto ref_seq = dag.GetReferenceSequence();
  auto id_to_cg_map = LoadCompactGenomeDataFromVCF(
      "data/test_ambiguous_vcf/SampleDAG_unique_ambiguous_leafs.vcf", ref_seq);

  SubtreeWeight<ParsimonyScore, MADAG> weight{dag};
  auto sample = AddMATConversion(weight.SampleTree({}));
  MAT::Tree mat;
  sample.View().BuildMAT(mat);

  // check these two methods apply the same set of changes
  MADAGApplyCompactGenomeData(dag_storage, id_to_cg_map);
  mat_apply_compact_genome_data(mat, id_to_cg_map);

  // this routine is from include/larch/impl/produce_mat.cpp
  check_MAT_MADAG_Eq(mat, dag);
}

[[maybe_unused]] void test_vcf_with_larch_usher() {
  std::string command =
      "./larch-usher -i data/test_ambiguous_vcf/unamb_mat.pb -r data / "
      "test_ambiguous_vcf / sample_reference_sequence.fasta - o "
      "test_larch_usher_output.pb -c 2 -v data / test_ambiguous_vcf / "
      "SampleDAG_unique_ambiguous_leafs.vcf ";

  assert_equal(0, std::system(command.c_str()), "Child process failed");
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[]() { test_ambiguous_vcf(); }, "Loading VCFs with Ambiguities Test"});

// [[maybe_unused]] static const auto test_added1 =
//     add_test({[]() { test_vcf_reading(); }, "load vcf and create MAT conversion"});

// [[maybe_unused]] static const auto test_added2 =
//     add_test({[]() {
//     test_vcf_compatible(); }, "Loading VCFs with Ambiguities on a
//     DAG and check pendant edge Mutations"});

// [[maybe_unused]] static const auto test_added3 =
//     add_test({
//   []() { test_vcf_with_larch_usher(); }, "Loading VCFs with Ambiguities
//                                              and running with MatOptimize "});

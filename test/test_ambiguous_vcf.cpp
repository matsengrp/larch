#include "test_common.hpp"
#include "sample_dag.hpp"

#include "larch/dag_loader.hpp"
#include "larch/spr/spr_view.hpp"
#include "larch/mat_conversion.hpp"
#include "larch/usher_glue.hpp"

std::ostream &operator<<(std::ostream &os, const MAT::Mutation mut) {
  os << mut.get_string();
  return os;
}

[[maybe_unused]] static std::string madag_info(const MADAGStorage<> &dag_storage) {
  std::stringstream os;
  auto dag = dag_storage.View();
  auto ref_seq = dag.GetReferenceSequence();
  os << "ref_seq: " << ref_seq << std::endl;
  os << "=== NODES: " << dag.GetNodesCount() << " ===" << std::endl;
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

[[maybe_unused]] static auto convert_madag_to_mat(const MADAGStorage<> &dag_storage) {
  MAT::Tree mat;
  auto dag = dag_storage.View();
  TestAssert(dag.IsTree());
  auto mat_conv = AddMATConversion(dag);
  mat_conv.View().BuildMAT(mat);
  return mat;
}

[[maybe_unused]] static auto load_mat_from_protobuf(std::string_view pb_path) {
  auto dag_storage = LoadDAGFromProtobuf(pb_path);
  auto mat = convert_madag_to_mat(dag_storage);
  return mat;
}

[[maybe_unused]] static void madag_label_nodes(MADAGStorage<> &dag_storage,
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

[[maybe_unused]] static std::string madag_to_fasta(const MADAGStorage<> &dag_storage,
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

[[maybe_unused]] static void madag_to_fasta(const MADAGStorage<> &dag_storage,
                                            const std::string &out_path, bool use_ids,
                                            bool leaves_only, bool include_reference) {
  std::ofstream file_out;
  file_out.open(out_path);
  file_out << madag_to_fasta(dag_storage, use_ids, leaves_only, include_reference)
           << std::endl;
  file_out.close();
}

[[maybe_unused]] static void fasta_to_vcf(const std::string &fasta_path_in,
                                          const std::string &vcf_path_out,
                                          bool include_reference = false) {
  std::string command, ref_command;
  command = "faToVcf --help";
  if (std::system(command.c_str()) != 0) {
    std::cerr << "ERROR: faToVcf installation required to call fasta_to_vcf."
              << std::endl;
    return;
  }
  ref_command = (include_reference ? "-ref=ref" : "");
  command = "faToVcf " + ref_command + " " + fasta_path_in + " " + vcf_path_out;
  if (std::system(command.c_str()) != 0) {
    std::cerr << "ERROR: faToVcf failed. Check your filepaths." << std::endl;
  }
}

[[maybe_unused]] static void mat_apply_compact_genome_data(
    MAT::Tree &tree, const std::unordered_map<std::string, CompactGenomeData> &mut_map,
    bool silence_warnings = true) {
  for (const auto &[name, muts] : mut_map) {
    auto node = tree.get_node(tree.node_name_to_node_idx(name));
    if (node) {
      node->mutations.clear();
      node->mutations.reserve(muts.size());
      for (const auto &[pos, base] : muts) {
        TestAssert(pos.value != NoId);
        MAT::Mutation mat_mut(
            "ref", static_cast<int>(pos.value), EncodeBaseMAT(base.ToChar()),
            EncodeBaseMAT(base.ToChar()), EncodeBaseMAT(base.ToChar()));
        node->mutations.push_back(mat_mut);
      }
    } else if (!silence_warnings) {
      std::cerr << "WARNING: Could not find sample_id `" << name << "` in MAT."
                << std::endl;
    }
  }
}

[[maybe_unused]] static bool madag_compare(const MADAGStorage<> &lhs_storage,
                                           const MADAGStorage<> &rhs_storage) {
  auto lhs = lhs_storage.View();
  auto rhs = rhs_storage.View();
  if (lhs.GetNodesCount() != rhs.GetNodesCount()) return false;
  if (lhs.GetEdgesCount() != rhs.GetEdgesCount()) return false;
  for (auto lhs_node : lhs.GetNodes()) {
    auto rhs_node = rhs.Get(lhs_node.GetId());
    if (lhs_node.GetCompactGenome() != rhs_node.GetCompactGenome()) return false;
    for (auto lhs_edge : lhs_node.GetChildren()) {
      auto rhs_edge = rhs.Get(lhs_edge.GetId());
      if (lhs_edge.GetEdgeMutations() != rhs_edge.GetEdgeMutations()) return false;
    }
  }
  return true;
}

bool operator==(const MADAGStorage<> &lhs_storage, const MADAGStorage<> &rhs_storage) {
  return madag_compare(lhs_storage, rhs_storage);
}

// MAIN TESTS

static bool silence_warnings = false;

[[maybe_unused]] static auto test_ambiguous_vcf() {
  // Reference Sequence (for when reading from pb)
  auto ref_seq = "GAA";
  // Fasta paths
  std::string amb_fasta_path = "data/test_ambiguous_vcf/amb.fasta";
  std::string unamb_fasta_path = "data/test_ambiguous_vcf/unamb.fasta";
  // VCF paths
  std::string amb_vcf_path = "data/test_ambiguous_vcf/amb.vcf";
  std::string unamb_vcf_path = "data/test_ambiguous_vcf/unamb.vcf";

  // Truth MADAG and MAT created from sample DAG.
  auto amb_dag_truth_storage = make_ambiguous_sample_dag();
  madag_label_nodes(amb_dag_truth_storage, false);
  auto amb_mat_truth = convert_madag_to_mat(amb_dag_truth_storage);
  auto unamb_dag_truth_storage = make_unambiguous_sample_dag();
  madag_label_nodes(unamb_dag_truth_storage, false);
  auto unamb_mat_truth = convert_madag_to_mat(unamb_dag_truth_storage);

  // Create fasta and vcf files from truth MADAGs.
  bool include_ref = true;
  madag_to_fasta(amb_dag_truth_storage, amb_fasta_path, false, true, include_ref);
  // fasta_to_vcf(amb_fasta_path, amb_vcf_path, include_ref);
  auto amb_cg = LoadCompactGenomeDataFromVCF(amb_vcf_path, ref_seq);
  madag_to_fasta(unamb_dag_truth_storage, unamb_fasta_path, false, true, include_ref);
  // fasta_to_vcf(unamb_fasta_path, unamb_vcf_path, include_ref);
  auto unamb_cg = LoadCompactGenomeDataFromVCF(unamb_vcf_path, ref_seq);

  // Create DAG with differing leaf sequences and apply VCFs that match truth DAGs.
  auto amb_dag_storage = make_unambiguous_sample_dag_2();
  madag_label_nodes(amb_dag_storage, false);
  TestAssert(not(amb_dag_storage == amb_dag_truth_storage) &&
             "Ambiguous DAGs found incorrectly equal before applying VCF.");
  LoadVCFData(amb_dag_storage, amb_vcf_path, silence_warnings);
  TestAssert((amb_dag_storage == amb_dag_truth_storage) &&
             "Ambiguous DAGs not equal after applying VCF.");
  auto unamb_dag_storage = make_unambiguous_sample_dag_2();
  madag_label_nodes(unamb_dag_storage, false);
  TestAssert(not(unamb_dag_storage == unamb_dag_truth_storage) &&
             "Unambiguous DAGs found incorrectly equal before applying VCF.");
  LoadVCFData(unamb_dag_storage, unamb_vcf_path, silence_warnings);
  TestAssert((unamb_dag_storage == unamb_dag_truth_storage) &&
             "Unambiguous DAGs not equal after applying VCF.");
}

void test_vcf_compatible() {
  auto dag_storage = make_sample_dag();
  madag_label_nodes(dag_storage, false);
  MADAG dag = dag_storage.View();
  auto ref_seq = dag.GetReferenceSequence();
  auto id_to_cg_map = LoadCompactGenomeDataFromVCF(
      "data/test_ambiguous_vcf/SampleDAG_unique_ambiguous_leafs.vcf", ref_seq);

  // make sure that each of the leaves in dag match exactly one key of id_to_cg_map
  for (auto leaf_id : dag.GetLeafs()) {
    if (leaf_id.GetSampleId().has_value()) {
      TestAssert(id_to_cg_map.find(leaf_id.GetSampleId().value()) !=
                 id_to_cg_map.end());
    }
  }

  // make sure that each parent edge is "compatible" with the leaf nodes
  MADAGApplyCompactGenomeData(dag_storage, id_to_cg_map, silence_warnings);
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
      TestAssert(leaf_cg.IsCompatible(parent_cg, ref_seq));
    }
  }
}

void test_vcf_reading() {
  auto dag_storage = make_sample_dag();
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
  MADAGApplyCompactGenomeData(dag_storage, id_to_cg_map, silence_warnings);
  mat_apply_compact_genome_data(mat, id_to_cg_map, silence_warnings);

  // this routine is from include/larch/impl/produce_mat.cpp
  check_MAT_MADAG_Eq(mat, dag);
}

[[maybe_unused]] void test_vcf_with_larch_usher() {
  std::string command =
      "./larch-usher -i data/test_ambiguous_vcf/unamb_mat.pb -r "
      "data/test_ambiguous_vcf/sample_reference_sequence.fasta -o "
      "test_larch_usher_output.pb -c 2 -v "
      "data/test_ambiguous_vcf/SampleDAG_unique_ambiguous_leafs.vcf ";

  TestAssert(0 == std::system(command.c_str()));
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[]() { test_ambiguous_vcf(); }, "Loading VCFs with Ambiguities Test"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[]() { test_vcf_reading(); }, "load vcf and create MAT conversion"});

[[maybe_unused]] static const auto test_added2 = add_test(
    {[]() { test_vcf_compatible(); },
     "Loading VCFs with Ambiguities on a DAG and check pendant edge Mutations "});

[[maybe_unused]] static const auto test_added3 =
    add_test({[]() { test_vcf_with_larch_usher(); },
              "Loading VCFs with Ambiguities and running with MatOptimize "});

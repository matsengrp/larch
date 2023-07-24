
#include "test_common.hpp"

#include "larch/dag_loader.hpp"
#include "larch/spr/spr_view.hpp"
#include "larch/mat_conversion.hpp"

#include "larch/usher_glue.hpp"

// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wold-style-cast"
// #pragma GCC diagnostic ignored "-Wunused-parameter"
// #pragma GCC diagnostic ignored "-Wdeprecated-copy-with-user-provided-copy"
// #pragma GCC diagnostic ignored "-Wunused-function"
// #pragma GCC diagnostic ignored "-Wshadow"
// #pragma GCC diagnostic ignored "-Wdeprecated-copy"
// #pragma GCC diagnostic ignored "-Wextra"
// #pragma GCC diagnostic ignored "-Wunused-function"
// #pragma GCC diagnostic ignored "-Wconversion"
// #pragma GCC diagnostic ignored "-Wmissing-field-initializers"
// #include "src/matOptimize/import_vcf_fast.cpp"
// #pragma GCC diagnostic pop

#include "larch/import_vcf.hpp"

[[maybe_unused]] static auto ConvertDAGToMAT(const MADAGStorage &dag_storage) {
  MAT::Tree mat;
  auto dag = dag_storage.View();
  Assert(dag.IsTree());
  auto mat_conv = AddMATConversion(dag);
  mat_conv.View().BuildMAT(mat);
  return mat;
}

[[maybe_unused]] static auto LoadMATFromProtobuf(std::string_view pb_path) {
  auto dag_storage = LoadDAGFromProtobuf(pb_path);
  auto mat = ConvertDAGToMAT(dag_storage);
  return mat;
}

std::ostream &operator<<(std::ostream &os, const MAT::Mutation mut) {
  os << mut.get_string();
  return os;
}

[[maybe_unused]] static std::string MADAGInfo(const MADAGStorage &dag_storage) {
  std::stringstream os;
  auto dag = dag_storage.View();
  auto ref_seq = dag.GetReferenceSequence();
  os << "ref_seq: " << ref_seq << std::endl;
  os << "== NODES ==" << std::endl;
  for (const auto node : dag.GetNodes()) {
    os << "NODE: " << node.GetId() << " " << node.GetSampleId().value_or("<no_name>")
       << std::endl;
    os << "children: [ ";
    for (const auto node_id : node.GetChildNodes()) {
      os << node_id << " ";
    }
    os << "] " << std::endl;
    os << "seq: " << node.GetCompactGenome().ToSequence(ref_seq) << std::endl;
    os << "muts: " << node.GetCompactGenome().ToString() << std::endl;
  }
  return os.str();
}

[[maybe_unused]] static std::string MATInfo(MAT::Tree &mat) {
  std::stringstream os;
  os << "newick: " << mat.get_newick_string(true, false, false, true) << std::endl;

  for (const auto node : mat.breadth_first_expansion()) {
    os << "node: " << node->node_id << " " << std::endl;
    auto node_name = (mat.get_node_name(node->node_id) != "")
                         ? mat.get_node_name(node->node_id)
                         : "<no_name>";
    os << "name: " << node_name << std::endl;
    os << "children: [ ";
    for (auto child : node->children) {
      os << child->node_id << " ";
    }
    os << "]" << std::endl;
    os << "muts: [ ";
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

[[maybe_unused]] static std::string OriginalStateInfo(
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

[[maybe_unused]] static void MADAGLabelNodes(MADAGStorage &dag_storage,
                                             bool leaves_only = true) {
  auto dag = dag_storage.View();
  size_t id = 0;
  for (auto node : dag.GetNodes()) {
    if (!leaves_only || node.IsLeaf()) {
      auto prefix = (node.IsLeaf() ? "x" : "y");
      prefix = (node.IsUA() ? "r" : prefix);
      node.SetSampleId({prefix + std::to_string(node.GetId().value)});
      id++;
    }
  }
}

[[maybe_unused]] static std::string MADAGToFasta(const MADAGStorage &dag_storage,
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

[[maybe_unused]] static void MADAGToFasta(const MADAGStorage &dag_storage,
                                          const std::string &out_path, bool use_ids,
                                          bool leaves_only, bool include_reference) {
  std::ofstream file_out;
  file_out.open(out_path);
  file_out << MADAGToFasta(dag_storage, use_ids, leaves_only, include_reference)
           << std::endl;
  file_out.close();
}

[[maybe_unused]] static std::vector<std::string> SplitString(const std::string &str,
                                                             char delimiter) {
  std::vector<std::string> substrings;
  std::string substring;

  for (char ch : str) {
    if (ch == delimiter) {
      substrings.push_back(substring);
      substring.clear();
    } else {
      substring += ch;
    }
  }

  substrings.push_back(substring);

  return substrings;
}

using CompactGenomeData = ContiguousMap<MutationPosition, MutationBase>;
[[maybe_unused]] static std::unordered_map<std::string, CompactGenomeData>
ReadVCFToCompactGenomeData(const std::string &path, const std::string &ref_seq) {
  std::unordered_map<std::string, CompactGenomeData> mut_map;

  std::ifstream in;
  in.open(path);
  if (!in) {
    std::cerr << "Failed to open file: " << path << std::endl;
    Assert(in);
  }

  std::map<std::string, size_t> name_col_map;
  std::map<size_t, std::string> col_name_map;
  std::string line;

  // parse header
  while (std::getline(in, line)) {
    // meta data
    if (line.substr(0, 2) == "##") {
      continue;
    }
    // column labels (capture node labels)
    size_t field_id = 0;
    std::string field;
    std::istringstream iss(line);
    if (line.substr(0, 1) == "#") {
      while (iss >> field) {
        field_id++;
        name_col_map[field] = field_id;
        col_name_map[field_id] = field;
        if (field_id <= 9) continue;
        mut_map[field] = CompactGenomeData();
      }
      break;
    }
  }

  const auto chrom_col_id = name_col_map["#CHROME"];
  const auto pos_col_id = name_col_map["POS"];
  const auto ref_col_id = name_col_map["REF"];
  const auto alt_col_id = name_col_map["ALT"];
  const auto fmt_col_id = name_col_map["FORMAT"];

  auto InsertMutation = [&mut_map, &ref_seq](const std::string &node_label,
                                             MutationPosition pos, MutationBase base) {
    if (base.ToChar() != ref_seq[pos.value]) {
      mut_map[node_label][pos] = base;
    }
  };

  // parse entries
  size_t entry_count = 0;
  while (std::getline(in, line)) {
    std::map<size_t, MutationBase> base_map;
    MutationPosition pos;
    MutationBase base;

    std::string chrom_name;
    size_t field_id = 0;
    std::string field;
    std::istringstream iss(line);
    while (iss >> field) {
      field_id++;
      if (field_id == chrom_col_id) {
        chrom_name = field;
        if (entry_count == 0) {
          mut_map[chrom_name] = CompactGenomeData();
        }
      }
      // Get mutation position
      if (field_id == pos_col_id) {
        pos = MutationPosition{size_t(std::stoi(field))};
      }
      // Get reference sequence base
      if (field_id == ref_col_id) {
        base = MutationBase(field[0]);
        base_map[0] = base;
        InsertMutation(chrom_name, pos, base);
      }
      // Get alternate mutations at position.
      if (field_id == alt_col_id) {
        size_t base_id = 1;
        for (const auto &base_char : SplitString(field, ',')) {
          base_map[base_id] = MutationBase(base_char[0]);
          base_id++;
        }
      }
      // Get samples from other nodes.
      if (field_id > fmt_col_id) {
        if (field == ".") {
          base = MutationBase('N');
          InsertMutation(col_name_map[field_id], pos, base);
        } else if (base_map.find(std::stoi(field)) != base_map.end()) {
          base = base_map[size_t(std::stoi(field))];
          InsertMutation(col_name_map[field_id], pos, base);
        } else {
          std::cerr << "Invalid symbol found on col " << col_name_map[field_id] << ": "
                    << field << std::endl;
          Assert(false);
        }
      }
    }
    entry_count++;
  }
  in.close();

  return mut_map;
}

[[maybe_unused]] static void MADAGApplyCompactGenomeData(
    MADAGStorage &dag_storage,
    const std::unordered_map<std::string, CompactGenomeData> &mut_map) {
  auto dag = dag_storage.View();
  // Convert node names to ids.
  std::unordered_map<NodeId, CompactGenomeData> tmp_mut_map;
  for (auto node : dag.GetNodes()) {
    if (node.GetSampleId().has_value() and
        mut_map.find(node.GetSampleId().value()) != mut_map.end()) {
      const auto &[name, muts] = *mut_map.find(node.GetSampleId().value());
      std::ignore = name;
      tmp_mut_map[node.GetId()] = CompactGenomeData();
      for (const auto &[pos, base] : muts) {
        tmp_mut_map[node.GetId()][pos] = base;
      }
    }
  }
  dag.SetCompactGenomesFromNodeMutationMap(std::move(tmp_mut_map));
  dag.RecomputeEdgeMutations();
}

[[maybe_unused]] static void MATApplyCompactGenomeData(
    MAT::Tree &tree,
    const std::unordered_map<std::string, CompactGenomeData> &mut_map) {
  for (const auto &[name, muts] : mut_map) {
    auto node = tree.get_node(name);
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

[[maybe_unused]] static void MATApplyVCF(MAT::Tree &tree, const std::string &vcf_path) {
  VCF_input(vcf_path.c_str(), tree);
}

[[maybe_unused]] static auto test_ambiguous_vcf() {
  // Reference Sequence
  auto ref_seq = "GAA";
  // Protobufs
  // std::string amb_pb = "data/test_ambiguous_vcf/amb.pb";
  // std::string unamb_pb = "data/test_ambiguous_vcf/unamb.pb";
  // VCFs
  std::string amb_vcf_path = "data/test_ambiguous_vcf/amb.vcf";
  std::string unamb_vcf_path = "data/test_ambiguous_vcf/unamb.vcf";

  // Create topology and

  // Build empty topology with labels (same topology as protobufs).
  // auto topo_dag_storage = MakeSampleDAGTopology();
  auto topo_dag_storage = MakeUnambiguousSampleDAG();
  MADAGLabelNodes(topo_dag_storage, false);

  // Make MAT from ambiguous sample DAG (same as amb.pb)
  auto amb_dag_storage = MakeAmbiguousSampleDAG();
  MADAGLabelNodes(amb_dag_storage, false);
  std::cout << "AMB_MADAG: " << std::endl;
  std::cout << MADAGInfo(amb_dag_storage) << std::endl;
  MADAGToFasta(amb_dag_storage, "data/test_ambiguous_vcf/amb.fasta", false, true,
               false);
  auto amb_mat_truth = ConvertDAGToMAT(amb_dag_storage);
  // Make MAT from unambiguous sample DAG (same as unamb.pb)
  auto unamb_dag_storage = MakeUnambiguousSampleDAG();
  MADAGLabelNodes(unamb_dag_storage, false);
  std::cout << "UNAMB_MADAG: " << std::endl;
  std::cout << MADAGInfo(unamb_dag_storage) << std::endl;
  MADAGToFasta(unamb_dag_storage, "data/test_ambiguous_vcf/unamb.fasta", false, true,
               false);
  auto unamb_mat_truth = ConvertDAGToMAT(unamb_dag_storage);

  // MAT from ambiguous VCF via MatOptimize.
  auto amb_mat_vcf = ConvertDAGToMAT(topo_dag_storage);
  MATApplyVCF(amb_mat_vcf, amb_vcf_path);
  // MAT from unambiguous VCF via MatOptimize.
  auto unamb_mat_vcf = ConvertDAGToMAT(topo_dag_storage);
  MATApplyVCF(unamb_mat_vcf, unamb_vcf_path);

  // MAT from ambiguous VCF via Compact Genome data.
  auto amb_mat_cg = ConvertDAGToMAT(topo_dag_storage);
  MATApplyCompactGenomeData(amb_mat_cg,
                            ReadVCFToCompactGenomeData(amb_vcf_path, ref_seq));
  // MAT from unambiguous VCF via Compact Genome data.
  auto unamb_mat_cg = ConvertDAGToMAT(topo_dag_storage);
  MATApplyCompactGenomeData(unamb_mat_cg,
                            ReadVCFToCompactGenomeData(unamb_vcf_path, ref_seq));

  // - create an Original_State_t object to hold the ambiguous leaf data
  Original_State_t amb_mat_truth_og, unamb_mat_truth_og;
  Original_State_t amb_mat_vcf_og, unamb_mat_vcf_og;
  Original_State_t amb_mat_cg_og, unamb_mat_cg_og;

  // - call check_samples() on the MAT and pass in the Original_State_t object to
  // populate it.
  check_samples(amb_mat_truth.root, amb_mat_truth_og, &amb_mat_truth, true);
  check_samples(unamb_mat_truth.root, unamb_mat_truth_og, &unamb_mat_truth, true);

  check_samples(amb_mat_vcf.root, amb_mat_vcf_og, &amb_mat_vcf, true);
  check_samples(unamb_mat_vcf.root, unamb_mat_vcf_og, &unamb_mat_vcf, true);

  check_samples(amb_mat_cg.root, amb_mat_cg_og, &amb_mat_cg, true);
  check_samples(unamb_mat_cg.root, unamb_mat_cg_og, &unamb_mat_cg, true);

  std::cout << "amb_mat_truth: " << OriginalStateInfo(amb_mat_truth_og) << std::endl;
  std::cout << "unamb_mat_truth: " << OriginalStateInfo(unamb_mat_truth_og)
            << std::endl;

  // std::cout << "amb_mat_vcf: " << OriginalStateInfo(amb_mat_vcf_og) << std::endl;
  // std::cout << "unamb_mat_vcf: " << OriginalStateInfo(unamb_mat_vcf_og) << std::endl;

  std::cout << "amb_mat_cg: " << OriginalStateInfo(amb_mat_cg_og) << std::endl;
  std::cout << "unamb_mat_cg: " << OriginalStateInfo(unamb_mat_cg_og) << std::endl;
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[]() { test_ambiguous_vcf(); }, "Loading VCFs with Ambiguities Test"});


#include "test_common.hpp"

#include "larch/dag_loader.hpp"
#include "larch/spr/spr_view.hpp"
#include "larch/mat_conversion.hpp"

#include "larch/usher_glue.hpp"

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
    os << "node: " << node.GetId() << std::endl;
    os << "name: " << node.GetSampleId().value_or("<no_name>") << std::endl;
    os << "children: [ ";
    for (const auto node_id : node.GetChildNodes()) {
      os << node_id << " ";
    }
    os << "] " << std::endl;
    os << "seq: " << node.GetCompactGenome().ToSequence(ref_seq) << std::endl;
    os << "muts: " << node.GetCompactGenome().ToString() << std::endl;
  }
  os << "== EDGES ==" << std::endl;
  for (const auto edge : dag.GetEdges()) {
    os << "edge: " << edge.GetId() << std::endl;
    os << "muts: " << edge.GetEdgeMutations().ToString() << std::endl;
  }
  return os.str();
}

[[maybe_unused]] static std::string MATTreeInfo(MAT::Tree &mat) {
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
ReadVCFToCompactGenomeData(const std::string &path) {
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

  const auto pos_col_id = name_col_map["POS"];
  const auto ref_col_id = name_col_map["REF"];
  const auto alt_col_id = name_col_map["ALT"];
  const auto fmt_col_id = name_col_map["FORMAT"];

  // parse entries
  while (std::getline(in, line)) {
    std::map<size_t, MutationBase> base_map;
    MutationPosition pos;
    MutationBase base;

    size_t field_id = 0;
    std::string field;
    std::istringstream iss(line);
    while (iss >> field) {
      field_id++;
      // Get mutation position
      if (field_id == pos_col_id) {
        pos = MutationPosition{size_t(std::stoi(field))};
      }
      // Get reference sequence base
      if (field_id == ref_col_id) {
        base_map[0] = MutationBase(field[0]);
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
          mut_map[col_name_map[field_id]][pos] = base;
        } else if (base_map.find(std::stoi(field)) != base_map.end()) {
          base = base_map[size_t(std::stoi(field))];
          mut_map[col_name_map[field_id]][pos] = base;
        } else {
          std::cerr << "Invalid symbol found on col " << col_name_map[field_id] << ": "
                    << field << std::endl;
          Assert(false);
        }
      }
    }
  }
  in.close();

  return mut_map;
}

[[maybe_unused]] static void MADAGApplyCompactGenome(
    MADAGStorage &dag_storage,
    const std::unordered_map<std::string, CompactGenomeData> &mut_map) {
  auto dag = dag_storage.View();
  // Convert node names to ids.
  std::unordered_map<NodeId, CompactGenomeData> tmp_mut_map;
  for (auto node : dag.GetNodes()) {
    if (mut_map.find(node.GetSampleId().value()) != mut_map.end()) {
      const auto &[name, muts] = *mut_map.find(node.GetSampleId().value());
      std::ignore = name;
      tmp_mut_map[node.GetId()] = CompactGenomeData();
      for (const auto &[pos, base] : muts) {
        tmp_mut_map[node.GetId()][pos] = base;
      }
    }
  }
  dag.SetCompactGenomesFromNodeMutationMap(std::move(tmp_mut_map));
}

[[maybe_unused]] static void MATApplyCompactGenome(
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

[[maybe_unused]] static auto test_ambiguous_vcf() {
  // std::string dag_pb_path = "data/test_ambiguous_vcf/unamb_dag.pb";
  // std::string mat_pb_path = "data/test_ambiguous_vcf/unamb_mat.pb";
  std::string amb_vcf_path = "data/test_ambiguous_vcf/amb.vcf";
  std::string unamb_vcf_path = "data/test_ambiguous_vcf/unamb.vcf";
  // std::string dag_pb_path2 = "data/test_ambiguous_vcf/unamb_dag2.pb";
  // std::string mat_pb_path2 = "data/test_ambiguous_vcf/unamb_mat2.pb";

  auto topo_dag_storage = MakeSampleDAGTopology();
  MADAGLabelNodes(topo_dag_storage, false);
  // std::cout << "TOPO_DAG: " << MADAGInfo(topo_dag_storage) << std::endl;
  // std::cout << "TOPO_FASTA: " << std::endl
  //           << MADAGToFasta(topo_dag_storage, false, false) << std::endl;

  auto amb_dag_storage = MakeAmbiguousSampleDAG();
  MADAGLabelNodes(amb_dag_storage, false);
  // MADAGToFasta(amb_dag_storage, "data/test_ambiguous_vcf/amb.fasta", false, true,
  //              false);

  auto unamb_dag_storage = MakeUnambiguousSampleDAG();
  MADAGLabelNodes(unamb_dag_storage, false);
  // MADAGToFasta(unamb_dag_storage, "data/test_ambiguous_vcf/unamb.fasta", false, true,
  //              false);

  auto mat_wo_vcf = ConvertDAGToMAT(topo_dag_storage);
  auto mat_w_vcf = ConvertDAGToMAT(topo_dag_storage);

  auto amb_cg_data = ReadVCFToCompactGenomeData(amb_vcf_path);
  auto unamb_cg_data = ReadVCFToCompactGenomeData(unamb_vcf_path);

  MATApplyCompactGenome(mat_w_vcf, unamb_cg_data);
  MATTreeInfo(mat_w_vcf);

  auto dag_w_vcf = MakeSampleDAGTopology();
  MADAGApplyCompactGenome(dag_w_vcf, unamb_cg_data);
  MADAGInfo(dag_w_vcf);

  // std::cout << "AMB_MAT_WO_VCF: " << MATTreeInfo(mat_wo_vcf) << std::endl;
  // MATApplyVCF(mat_w_vcf, unamb_vcf_path);
  // std::cout << "AMB_MAT_W_VCF: " << MATTreeInfo(mat_w_vcf) << std::endl;

  // Original_State_t amb_og_state;
  // check_samples(amb_mat_w_vcf.root, amb_og_state, &amb_mat_w_vcf, true);
  // std::cout << "OG_STATE: " << OriginalStateInfo(amb_og_state) << std::endl;

  // auto unamb_dag_storage = MakeUnambiguousSampleDAG();
  // auto unamb_dag = unamb_dag_storage.View();
  // Assert(unamb_dag.IsTree());
  // // std::cout << "UNAMB_DAG: " << MADAGInfo(unamb_dag_storage) << std::endl;
  // auto unamb_mat = ConvertDAGToMAT(unamb_dag_storage);
  // std::cout << "UNAMB_MAT: " << MATTreeInfo(unamb_mat) << std::endl;

  // // - construct a DAG from protobuf (tree shaped, since it's from a MAT protobuf)
  // auto dag_storage = LoadDAGFromProtobuf(dag_pb_path);
  // auto dag = dag_storage.View();
  // Assert(dag.IsTree());
  // // MADAGLabelNodes(dag_storage);
  // dag.RecomputeCompactGenomes(true);
  // std::cout << "DAG_FROM_PB: " << MADAGInfo(dag_storage) << std::endl;
  // std::cout << "FASTA" << std::endl << MADAGToFastaString(dag_storage) <<
  // std::endl;
  // // StoreDAGToProtobuf(dag_storage.View(), dag_pb_path2);
  // // StoreTreeToProtobuf(dag_storage.View(), mat_pb_path2);

  // // auto tree_storage =
  // //     LoadTreeFromProtobuf(mat_pb_path,
  // dag_storage.View().GetReferenceSequence());
  // // auto tree = tree_storage.View();
  // // Assert(tree.IsTree());
  // // tree.RecomputeCompactGenomes(true);
  // // std::cout << "TREE_FROM_PB: " << MADAGInfo(tree_storage) << std::endl;

  // auto mat_from_pb = LoadMATFromProtobuf(dag_pb_path);
  // auto mat_from_vcf = LoadMATFromProtobuf(dag_pb_path);

  // // - call VCF_input() on the MAT
  // ApplyVCFToMAT(mat_from_vcf, vcf_path);
  // std::cout << "MAT_FROM_PB: " << MATTreeInfo(mat_from_pb) << std::endl;
  // std::cout << "MAT_FROM_VCF: " << MATTreeInfo(mat_from_vcf) << std::endl;

  // - create an Original_State_t object to hold the ambiguous leaf data
  Original_State_t og_state;
  // - call check_samples() on the MAT and pass in the Original_State_t object to
  // populate it.
  check_samples(mat_w_vcf.root, og_state, &mat_w_vcf, true);
  std::cout << "og_state [after]: " << OriginalStateInfo(og_state) << std::endl;
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[]() { test_ambiguous_vcf(); }, "Loading VCFs with Ambiguities Test"});

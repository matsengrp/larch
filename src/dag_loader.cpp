#include <fstream>
#include <vector>
#include <unordered_map>
#include <iostream>

#include <sys/stat.h>
#include <fcntl.h>

#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wredundant-decls"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Woverflow"
#pragma GCC diagnostic ignored "-Wpedantic"
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/gzip_stream.h>

#pragma GCC diagnostic ignored "-Wstack-usage="
#include "nlohmann/json.hpp"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>

#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/dag_loader.hpp"
#include "proto/dag.pb.h"
#include "parsimony.pb.h"
#include "larch/newick.hpp"
#include "larch/dagbin_fileio.hpp"

namespace {

auto OpenFile(std::string_view path) {
  std::ifstream in{std::string{path}};
  if (!in) {
    std::cerr << "Failed to open file: '" << path << "'." << std::endl;
    Assert(in);
  }
  return in;
}

bool IsGzipped(std::string_view path) {
  std::ifstream in = OpenFile(path);
  std::array<unsigned char, 2> header{};
  in >> header[0] >> header[1];
  constexpr const unsigned char HeaderMagic0 = 0x1f;
  constexpr const unsigned char HeaderMagic1 = 0x8b;
  return header[0] == HeaderMagic0 and header[1] == HeaderMagic1;
}

template <typename T>
void Parse(T& data, std::string_view path) {
  if (IsGzipped(path)) {
    google::protobuf::io::FileInputStream in_compressed{
        open(std::string{path}.c_str(), O_RDONLY)};
    in_compressed.SetCloseOnDelete(true);
    google::protobuf::io::GzipInputStream in{&in_compressed};
    [[maybe_unused]] bool parsed = data.ParseFromZeroCopyStream(&in);
    Assert(parsed);
  } else {
    std::ifstream in = OpenFile(path);
    [[maybe_unused]] bool parsed = data.ParseFromIstream(&in);
    Assert(parsed);
  }
}

}  // namespace

MADAGStorage<> LoadDAG(std::string_view input_dag_path, FileFormat file_format,
                       std::optional<std::string> refseq_path) {
  if (file_format == FileFormat::Infer) {
    file_format = InferFileFormat(input_dag_path);
  }
  auto load_dag = [&]() {
    switch (file_format) {
      case FileFormat::Dagbin:
        return LoadDAGFromDagbin(input_dag_path);
      case FileFormat::ProtobufDAG:
        return LoadDAGFromProtobuf(input_dag_path);
      case FileFormat::ProtobufTree:
        return LoadTreeFromProtobuf(
            input_dag_path,
            refseq_path.has_value() ? LoadReferenceSequence(refseq_path.value()) : "");
      case FileFormat::JsonDAG:
        return LoadDAGFromJson(input_dag_path);
      case FileFormat::Protobuf:
      case FileFormat::DebugAll:
      default:
        std::cerr
            << "ERROR: Could not load DAG with unrecognized/unsupported file format '"
            << input_dag_path << "'." << std::endl;
        std::exit(EXIT_FAILURE);
    }
  };
  return load_dag();
}

MADAGStorage<> LoadDAGFromDagbin(std::string_view path) {
  return DagbinFileIO::ReadDAG(path);
}

MADAGStorage<> LoadDAGFromProtobuf(std::string_view path) {
  ProtoDAG::data data;
  Parse(data, path);

  MADAGStorage<> result_storage = MADAGStorage<>::EmptyDefault();
  auto result = result_storage.View();
  result.SetReferenceSequence(data.reference_seq());

  for (const auto& i : data.node_names()) {
    auto new_node = result.AddNode({static_cast<size_t>(i.node_id())});
    if (i.condensed_leaves().size() > 0) {
      for (auto cl : i.condensed_leaves()) {
        new_node = SampleId{cl};
        break;
      }
    }
  }

  size_t edge_id = 0;
  for (const auto& i : data.edges()) {
    auto edge = result.AddEdge({edge_id++}, {static_cast<size_t>(i.parent_node())},
                               {static_cast<size_t>(i.child_node())},
                               {static_cast<size_t>(i.parent_clade())});
    EdgeMutations muts;
    for (const auto& mut : i.edge_mutations()) {
      static const std::array<char, 4> decode = {'A', 'C', 'G', 'T'};
      Assert(mut.position() > 0);
      Assert(mut.mut_nuc().size() == 1);
      muts[{static_cast<size_t>(mut.position())}] = {
          decode.at(static_cast<size_t>(mut.par_nuc())),
          decode.at(static_cast<size_t>(mut.mut_nuc().Get(0)))};
    }
    edge.SetEdgeMutations(std::move(muts));
  }
  result.BuildConnections();
  result.AssertUA();
  // result.RecomputeCompactGenomes();

  for (auto node : result.GetNodes()) {
    if (node.IsLeaf() and not node.HaveSampleId()) {
      node = SampleId{node.GetCompactGenome().ToString()};
    }
  }

  result.GetRoot().Validate(true, true);
  return result_storage;
}

// NOLINTNEXTLINE (cppcoreguidelines-interfaces-global-init)
static const auto DecodeMutation =
    [](auto& mut) -> std::pair<MutationPosition, std::pair<char, char>> {
  static const std::array<char, 4> decode = {'A', 'C', 'G', 'T'};
  Assert(mut.mut_nuc().size() == 1);
  return {{static_cast<size_t>(mut.position())},
          {decode.at(static_cast<size_t>(mut.par_nuc())),
           decode.at(static_cast<size_t>(mut.mut_nuc().Get(0)))}};
};

MADAGStorage<> LoadTreeFromProtobuf(std::string_view path,
                                    std::string_view reference_sequence) {
  Parsimony::data data;
  Parse(data, path);

  MADAGStorage<> result_storage = MADAGStorage<>::EmptyDefault();
  auto result = result_storage.View();
  result.SetReferenceSequence(reference_sequence);

  std::unordered_map<size_t, size_t> num_children;
  std::map<size_t, std::optional<std::string>> seq_ids;
  ParseNewick(
      data.newick(),
      [&seq_ids](size_t node_id, std::string_view label, std::optional<double>) {
        seq_ids[node_id] = label;
      },
      [&result, &num_children](size_t parent, size_t child) {
        result.AddEdge({child}, {parent}, {child}, {num_children[parent]++});
      });
  result.InitializeNodes(result.GetEdgesCount() + 1);
  result.BuildConnections();

  for (auto node : result.GetNodes()) {
    if (node.IsLeaf()) {
      node = SampleId{seq_ids[node.GetId().value]};
    }
  }

  result.AddUA({});

  Assert(static_cast<size_t>(data.node_mutations_size()) == result.GetNodesCount() - 1);
  auto apply_mutations = [](auto& self, auto edge, const auto& node_mutations,
                            size_t& idx) -> void {
    const auto& pb_muts = node_mutations.Get(static_cast<int>(idx++)).mutation();
    EdgeMutations muts;
    for (auto i : pb_muts | ranges::views::transform(DecodeMutation)) {
      muts.insert(i);
    }
    edge.SetEdgeMutations(std::move(muts));
    for (auto child : edge.GetChild().GetChildren()) {
      self(self, child, node_mutations, idx);
    }
  };
  size_t muts_idx = 0;
  apply_mutations(apply_mutations, result.GetRoot().GetFirstChild(),
                  data.node_mutations(), muts_idx);

  // uncollapsing has to occur _after_ we read the mutations, since those mutations
  // are stored for edges in an ordered traversal of the condensed newick tree.
  const auto& leaf_ids = result.GetNodes() | Transform::GetId();
  for (auto orig_leaf_id : leaf_ids) {
    if (result.Get(orig_leaf_id).IsLeaf()) {
      auto orig_leaf_in_dag = result.Get(orig_leaf_id);
      std::string node_name = orig_leaf_in_dag.GetSampleId().value();
      // check if this leaf is condensed
      if (node_name.find("_condensed_") != std::string::npos) {
        // find the corresponding condensed leaf in the data.condensed_nodes()
        // dictionary
        for (const auto& cn : data.condensed_nodes()) {
          std::string condensed_node_name = static_cast<std::string>(cn.node_name());
          if (condensed_node_name == node_name) {
            // expand this condensed leaf node in the result by:
            // - renaming this node to the first string in the list, and
            // - adding all of the sibling nodes as children of the condensed leaf
            // node's parent.
            auto parent_edge = orig_leaf_in_dag.GetSingleParent();
            auto parent_node = parent_edge.GetParent();
            auto clade_idx = parent_node.GetCladesCount();
            size_t ctr = 0;
            for (const auto& sib_node : cn.condensed_leaves()) {
              auto sib_node_name = static_cast<std::string>(sib_node);
              if (ctr < 1) {
                orig_leaf_in_dag = SampleId{sib_node_name};
              } else {
                auto muts_copy = parent_edge.GetEdgeMutations().Copy(&parent_edge);
                auto new_sib_node = result.AppendNode();
                new_sib_node = SampleId{sib_node_name};
                auto new_sib_edge = result.AppendEdge();
                new_sib_edge.SetEdgeMutations(std::move(muts_copy));
                result.AddEdge(new_sib_edge, parent_node, new_sib_node, {clade_idx++});
                result.AddLeaf(new_sib_node);
              }
              ctr++;
            }
          }
        }
      }
    }
  }

  result.BuildConnections();
  result.GetRoot().Validate(true, false);
  return result_storage;
}

[[nodiscard]] nlohmann::json LoadJson(std::string_view path) {
  if (IsGzipped(path)) {
    google::protobuf::io::FileInputStream in_compressed{
        open(std::string{path}.c_str(), O_RDONLY)};
    in_compressed.SetCloseOnDelete(true);
    google::protobuf::io::GzipInputStream in{&in_compressed};
    std::vector<char> bytes;
    const void* data{};
    int size{};
    while (in.Next(&data, &size)) {
      bytes.insert(bytes.end(), static_cast<const char*>(data),
                   static_cast<const char*>(data) + size);
    }
    return nlohmann::json::parse(bytes);
  }

  nlohmann::json result;
  std::ifstream in = OpenFile(path);
  in >> result;
  return result;
}

namespace {

CompactGenome GetCompactGenome(const nlohmann::json& json,
                               size_t compact_genome_index) {
  ContiguousMap<MutationPosition, MutationBase> result;
  result.reserve(json["compact_genomes"][compact_genome_index].size());
  for (const auto& mutation : json["compact_genomes"][compact_genome_index]) {
    MutationPosition position = {mutation[0]};
    std::string mut_nuc = mutation[1][1].get<std::string>();
    Assert(mut_nuc.size() == 1);
    result.insert({position, mut_nuc.at(0)});
  }
  return CompactGenome{std::move(result)};
}

}  // namespace

/*

compact_genome_list is a sorted list of compact genomes, where each compact
genome is a sorted list of mutations. Each mutation has the format
[seq_idx, [old_base, new_base]], where seq_idx is the 1-indexed nucleotide
sequence site of the mutation.

node_list is a list of [label_idx, clade_list] pairs, where label_idx is the
index of the node's compact genome in compact_genome_list, and clade_list is a
list of lists of compact_genome_list indices, encoding sets of child clades.

edge_list is a list of triples [parent_idx, child_idx, clade_idx], where
parent_idx is the index of the edge's parent node in node_list, child_idx is
the index of the edge's child node in node_list, and clade_idx is the index of
the clade in the parent node's clade_list from which this edge descends.

*/

MADAGStorage<> LoadDAGFromJson(std::string_view path) {
  nlohmann::json json = LoadJson(path);
  MADAGStorage<> result_storage = MADAGStorage<>::EmptyDefault();
  auto result = result_storage.View();
  result.SetReferenceSequence(std::string(json["refseq"][1]));

  size_t id = 0;
  for ([[maybe_unused]] auto& i : json["nodes"]) {
    auto node = result.AddNode({id++});
    size_t compact_genome_index = i[0];
    node = GetCompactGenome(json, compact_genome_index);
  }
  id = 0;
  for (auto& i : json["edges"]) {
    result.AddEdge({id++}, {i[0]}, {i[1]}, {i[2]});
  }
  result.BuildConnections();
  result.AssertUA();

  for (auto node : result.GetNodes()) {
    if (node.IsLeaf()) {
      node = SampleId{node.GetCompactGenome().ToString()};
    }
  }

  result.GetRoot().Validate(true, true);
  return result_storage;
}

std::string LoadReferenceSequence(std::string_view path) {
  std::string result;
  if (IsGzipped(path)) {
    std::ifstream file{std::string{path}};
    boost::iostreams::filtering_ostream str;
    str.push(boost::iostreams::gzip_decompressor());
    str.push(boost::iostreams::back_inserter(result));
    boost::iostreams::copy(file, str);
  } else {
    std::ifstream file{std::string{path}};
    while (file >> result) {
    }
  }
  return result;
}

namespace {

template <typename Mutation>
void InitMutation(Mutation* proto_mut, size_t pos, char ref, char par, char mut) {
  proto_mut->set_position(static_cast<int32_t>(pos));
  proto_mut->set_ref_nuc(EncodeBasePB(ref));
  proto_mut->set_par_nuc(EncodeBasePB(par));
  proto_mut->add_mut_nuc(EncodeBasePB(mut));
}

std::string ToEdgeMutationsString(const MAT::Node* node) {
  static const std::array<char, 4> decode = {'A', 'C', 'G', 'T'};
  std::string result = "<";
  for (const MAT::Mutation& mut : node->mutations) {
    result += decode.at(one_hot_to_two_bit(mut.get_par_one_hot()));
    result += std::to_string(mut.get_position());
    result += decode.at(one_hot_to_two_bit(mut.get_mut_one_hot()));
    result += ", ";
  }
  return result + ">";
}

template <typename iostream>
void MATToDOT(const MAT::Node* node, iostream& out,
              std::set<const MAT::Node*> visited) {
  Assert(visited.insert(node).second);

  for (auto* i : node->children) {
    MATToDOT(i, out, visited);
    out << "  \"" << node->node_id << " " << ToEdgeMutationsString(node) << "\" -> \""
        << i->node_id << " " << ToEdgeMutationsString(i) << "\"";
    // out << "[ headlabel=\"";
    // out << ToEdgeMutationsString(i);
    // out << "\" ]";
    out << "\n";
  }
}

}  // namespace

template <typename iostream>
void MATToDOT(const MAT::Tree& mat, iostream& out) {
  out << "digraph {\n";
  out << "  forcelabels=true\n";
  out << "  nodesep=1.0\n";
  out << "  ranksep=2.0\n";
  out << "  ratio=1.0\n";
  std::set<const MAT::Node*> visited;
  MATToDOT(mat.root, out, visited);
  out << "}\n";
}

std::unordered_map<std::string, CompactGenomeData> LoadCompactGenomeDataFromVCF(
    const std::string& path, const std::string& ref_seq) {
  std::unordered_map<std::string, CompactGenomeData> mut_map;

  std::ifstream in = OpenFile(path);

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

  const auto chrom_col_id = name_col_map["#CHROM"];
  const auto pos_col_id = name_col_map["POS"];
  const auto ref_col_id = name_col_map["REF"];
  const auto alt_col_id = name_col_map["ALT"];
  const auto fmt_col_id = name_col_map["FORMAT"];

  auto InsertMutation = [&mut_map, &ref_seq](const std::string& node_label,
                                             MutationPosition pos, MutationBase base) {
    if (base.ToChar() != ref_seq[pos.value - 1]) {
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
      // Get chrom name
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
        for (const auto& base_char : SplitString(field, ',')) {
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

void MADAGApplyCompactGenomeData(
    MADAGStorage<>& dag_storage,
    const std::unordered_map<std::string, CompactGenomeData>& mut_map,
    bool silence_warnings) {
  auto dag = dag_storage.View();
  // Convert node names to ids.
  std::unordered_map<std::string, bool> visited_ids;
  NodeMutMap tmp_mut_map;
  for (auto node : dag.GetNodes()) {
    if (node.GetSampleId().has_value() and
        mut_map.find(node.GetSampleId().value()) != mut_map.end()) {
      const auto& [name, muts] = *mut_map.find(node.GetSampleId().value());
      if (!silence_warnings) {
        visited_ids[name] = true;
      }
      tmp_mut_map[node.GetId()] = CompactGenomeData();
      for (const auto& [pos, base] : muts) {
        tmp_mut_map[node.GetId()][pos] = base;
      }
    }
  }
  if (!silence_warnings) {
    for (auto& name_muts : mut_map) {
      if (visited_ids.find(name_muts.first) == visited_ids.end()) {
        std::cerr << "WARNING: Could not find sample_id `" << name_muts.first
                  << "` in MADAG." << std::endl;
      }
    }
  }
  dag.UpdateCompactGenomesFromNodeMutationMap(std::move(tmp_mut_map));
  dag.RecomputeEdgeMutations();
}

void LoadVCFData(MADAGStorage<>& dag_storage, std::string& vcf_path,
                 bool silence_warnings) {
  if (not vcf_path.empty()) {
    auto ref_seq = dag_storage.View().GetReferenceSequence();
    auto cg_data = LoadCompactGenomeDataFromVCF(vcf_path, ref_seq);
    MADAGApplyCompactGenomeData(dag_storage, cg_data, silence_warnings);
  }
}

#include <fstream>
#include <vector>
#include <unordered_map>
#include <ostream>

#include <sys/stat.h>
#include <fcntl.h>

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/gzip_stream.h>

#pragma GCC diagnostic ignored "-Wstack-usage="
#include "nlohmann/json.hpp"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>

#include "larch/dag_loader.hpp"
#include "dag.pb.h"
#include "parsimony.pb.h"
#include "larch/newick.hpp"

static bool IsGzipped(std::string_view path) {
  std::ifstream in{std::string{path}};
  Assert(in);
  std::array<unsigned char, 2> header{};
  in >> header[0] >> header[1];
  constexpr const unsigned char HeaderMagic0 = 0x1f;
  constexpr const unsigned char HeaderMagic1 = 0x8b;
  return header[0] == HeaderMagic0 and header[1] == HeaderMagic1;
}

template <typename T>
static void Parse(T& data, std::string_view path) {
  if (IsGzipped(path)) {
    google::protobuf::io::FileInputStream in_compressed{
        open(std::string{path}.c_str(), O_RDONLY)};
    in_compressed.SetCloseOnDelete(true);
    google::protobuf::io::GzipInputStream in{&in_compressed};
    bool parsed = data.ParseFromZeroCopyStream(&in);
    Assert(parsed);
  } else {
    std::ifstream in{std::string{path}};
    bool parsed = data.ParseFromIstream(&in);
    Assert(parsed);
  }
}

MADAGStorage LoadDAGFromProtobuf(std::string_view path) {
  ProtoDAG::data data;
  Parse(data, path);

  MADAGStorage result;
  result.View().SetReferenceSequence(data.reference_seq());

  for (const auto& i : data.node_names()) {
    result.View().AddNode({static_cast<size_t>(i.node_id())});
  }

  size_t edge_id = 0;
  for (const auto& i : data.edges()) {
    auto edge = result.View().AddEdge(
        {edge_id++}, {static_cast<size_t>(i.parent_node())},
        {static_cast<size_t>(i.child_node())}, {static_cast<size_t>(i.parent_clade())});
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
  result.View().BuildConnections();
  result.View().AssertUA();
  return result;
}

static const auto DecodeMutation =
    [](auto& mut) -> std::pair<MutationPosition, std::pair<char, char>> {
  static const std::array<char, 4> decode = {'A', 'C', 'G', 'T'};
  Assert(mut.mut_nuc().size() == 1);
  return {{static_cast<size_t>(mut.position())},
          {decode.at(static_cast<size_t>(mut.par_nuc())),
           decode.at(static_cast<size_t>(mut.mut_nuc().Get(0)))}};
};

MADAGStorage LoadTreeFromProtobuf(std::string_view path,
                                  std::string_view reference_sequence) {
  Parsimony::data data;
  Parse(data, path);

  MADAGStorage result;
  result.View().SetReferenceSequence(reference_sequence);

  std::unordered_map<size_t, size_t> num_children;
  std::map<size_t, std::optional<std::string>> seq_ids;
  ParseNewick(
      data.newick(),
      [&seq_ids](size_t node_id, std::string_view label, std::optional<double>) {
        seq_ids[node_id] = label;
      },
      [&result, &num_children](size_t parent, size_t child) {
        result.View().AddEdge({child}, {parent}, {child}, {num_children[parent]++});
      });
  result.View().InitializeNodes(result.View().GetEdgesCount() + 1);
  result.View().BuildConnections();

  for (auto node : result.View().GetNodes()) {
    if (node.IsLeaf()) {
      node.SetSampleId(seq_ids[node.GetId().value]);
    }
  }

  result.View().AddUA({});

  Assert(static_cast<size_t>(data.node_mutations_size()) ==
         result.View().GetNodesCount() - 1);
  using Edge = MutableMADAG::EdgeView;
  auto apply_mutations = [](auto& self, Edge edge, const auto& node_mutations,
                            size_t& idx) -> void {
    const auto& pb_muts = node_mutations.Get(static_cast<int>(idx++)).mutation();
    EdgeMutations muts;
    for (auto i : pb_muts | ranges::views::transform(DecodeMutation)) {
      muts.insert(i);
    }
    edge.SetEdgeMutations(std::move(muts));
    for (Edge child : edge.GetChild().GetChildren()) {
      self(self, child, node_mutations, idx);
    }
  };
  size_t muts_idx = 0;
  apply_mutations(apply_mutations, result.View().GetRoot().GetFirstChild(),
                  data.node_mutations(), muts_idx);

  return result;
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
  std::ifstream in{std::string{path}};
  Assert(in);
  in >> result;
  return result;
}

static CompactGenome GetCompactGenome(const nlohmann::json& json,
                                      size_t compact_genome_index) {
  ContiguousMap<MutationPosition, char> result;
  result.reserve(json["compact_genomes"][compact_genome_index].size());
  for (const auto& mutation : json["compact_genomes"][compact_genome_index]) {
    MutationPosition position = {mutation[0]};
    std::string mut_nuc = mutation[1][1].get<std::string>();
    Assert(mut_nuc.size() == 1);
    result.insert({position, mut_nuc.at(0)});
  }
  return result;
}

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

MADAGStorage LoadDAGFromJson(std::string_view path) {
  nlohmann::json json = LoadJson(path);
  MADAGStorage result;
  result.View().SetReferenceSequence(std::string(json["refseq"][1]));

  size_t id = 0;
  for ([[maybe_unused]] auto& i : json["nodes"]) {
    auto node = result.View().AddNode({id++});
    size_t compact_genome_index = i[0];
    node = GetCompactGenome(json, compact_genome_index);
  }
  id = 0;
  for (auto& i : json["edges"]) {
    result.View().AddEdge({id++}, {i[0]}, {i[1]}, {i[2]});
  }
  result.View().BuildConnections();
  result.View().AssertUA();
  return result;
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

static void MATToDOT(const MAT::Node* node, std::ostream& out,
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

void MATToDOT(const MAT::Tree& mat, std::ostream& out) {
  out << "digraph {\n";
  out << "  forcelabels=true\n";
  out << "  nodesep=1.0\n";
  out << "  ranksep=2.0\n";
  out << "  ratio=1.0\n";
  std::set<const MAT::Node*> visited;
  MATToDOT(mat.root, out, visited);
  out << "}\n";
}

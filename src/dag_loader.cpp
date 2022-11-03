#include <fstream>
#include <vector>
#include <unordered_map>

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
  unsigned char header[2];
  in >> header[0] >> header[1];
  return header[0] == 0x1f and header[1] == 0x8b;
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

  for (auto& i : data.node_names()) {
    result.View().AddNode({static_cast<size_t>(i.node_id())});
  }

  size_t edge_id = 0;
  for (auto& i : data.edges()) {
    auto edge = result.View().AddEdge(
        {edge_id++}, {static_cast<size_t>(i.parent_node())},
        {static_cast<size_t>(i.child_node())}, {static_cast<size_t>(i.parent_clade())});
    auto& muts = edge.GetEdgeMutations();
    for (auto& mut : i.edge_mutations()) {
      static const char decode[] = {'A', 'C', 'G', 'T'};
      Assert(mut.position() > 0);
      Assert(mut.mut_nuc().size() == 1);
      muts[{static_cast<size_t>(mut.position())}] = {decode[mut.par_nuc()],
                                                     decode[mut.mut_nuc().Get(0)]};
    }
  }
  result.View().BuildConnections();
  result.View().AssertUA();
  return result;
}

static auto DecodeMutation =
    [](auto& mut) -> std::pair<MutationPosition, std::pair<char, char>> {
  static const char decode[] = {'A', 'C', 'G', 'T'};
  Assert(mut.mut_nuc().size() == 1);
  return {{static_cast<size_t>(mut.position())},
          {decode[mut.par_nuc()], decode[mut.mut_nuc().Get(0)]}};
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
      [&seq_ids](size_t node_id, std::string label, std::optional<double>) {
        seq_ids[node_id] = label;
      },
      [&result, &num_children](size_t parent, size_t child) {
        result.View().AddEdge({child}, {parent}, {child}, {num_children[parent]++});
      });
  result.View().InitializeNodes(result.View().GetEdgesCount() + 1);
  result.View().BuildConnections();

  for (auto node : result.View().GetNodes()) {
    if (node.IsLeaf()) {
      node.SetSampleId(std::move(seq_ids[node.GetId().value]));
    }
  }

  result.View().AddUA({});

  Assert(static_cast<size_t>(data.node_mutations_size()) ==
         result.View().GetNodesCount() - 1);

  auto apply_mutations = [](auto& self, MutableMADAG::Edge edge,
                            const auto& node_mutations, size_t& idx) -> void {
    const auto& pb_muts = node_mutations.Get(static_cast<int>(idx++)).mutation();
    for (auto i : pb_muts | ranges::views::transform(DecodeMutation)) {
      edge.GetEdgeMutations().insert(i);
    }
    for (MutableMADAG::Edge child : edge.GetChild().GetChildren()) {
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
    const void* data;
    int size;
    while (in.Next(&data, &size)) {
      bytes.insert(bytes.end(), static_cast<const char*>(data),
                   static_cast<const char*>(data) + size);
    }
    return nlohmann::json::parse(bytes);
  } else {
    nlohmann::json result;
    std::ifstream in{std::string{path}};
    Assert(in);
    in >> result;
    return result;
  }
}

static CompactGenome GetCompactGenome(const nlohmann::json& json,
                                      size_t compact_genome_index) {
  std::vector<std::pair<MutationPosition, char>> result;
  result.reserve(json["compact_genomes"][compact_genome_index].size());
  for (auto& mutation : json["compact_genomes"][compact_genome_index]) {
    MutationPosition position = {mutation[0]};
    std::string mut_nuc = mutation[1][1].get<std::string>();
    Assert(mut_nuc.size() == 1);
    result.emplace_back(position, mut_nuc.at(0));
  }
  std::sort(result.begin(), result.end(),
            [](auto lhs, auto rhs) { return lhs.first < rhs.first; });
  result.erase(std::unique(result.begin(), result.end(),
                           [](auto lhs, auto rhs) { return lhs.first == rhs.first; }),
               result.end());
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
    MutableMADAG::Node node = result.View().AddNode({id++});
    size_t compact_genome_index = i[0];
    node.SetCompactGenome(GetCompactGenome(json, compact_genome_index));
  }
  id = 0;
  for (auto& i : json["edges"]) {
    result.View().AddEdge({id++}, {i[0]}, {i[1]}, {i[2]});
  }
  result.View().BuildConnections();
  result.View().AssertUA();
  return result;
}

static int32_t EncodeBase(char base) {
  switch (base) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
    default:
      Fail("Invalid base");
  };
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
  proto_mut->set_ref_nuc(EncodeBase(ref));
  proto_mut->set_par_nuc(EncodeBase(par));
  proto_mut->add_mut_nuc(EncodeBase(mut));
}

void StoreDAGToProtobuf(MADAG dag, std::string_view path) {
  dag.AssertUA();
  ProtoDAG::data data;

  data.set_reference_seq(dag.GetReferenceSequence());

  for (size_t i = 0; i < dag.GetNodesCount(); ++i) {
    auto* proto_node = data.add_node_names();
    proto_node->set_node_id(static_cast<int64_t>(i));
  }

  for (MADAG::Edge edge : dag.GetEdges()) {
    auto* proto_edge = data.add_edges();
    proto_edge->set_edge_id(static_cast<int64_t>(edge.GetId().value));
    proto_edge->set_parent_node(static_cast<int64_t>(edge.GetParentId().value));
    proto_edge->set_child_node(static_cast<int64_t>(edge.GetChildId().value));
    proto_edge->set_parent_clade(static_cast<int64_t>(edge.GetClade().value));
    for (auto [pos, nucs] : edge.GetEdgeMutations()) {
      auto* proto_mut = proto_edge->add_edge_mutations();
      proto_mut->set_position(static_cast<int32_t>(pos.value));
      proto_mut->set_par_nuc(EncodeBase(nucs.first));
      proto_mut->add_mut_nuc(EncodeBase(nucs.second));
    }
  }

  std::ofstream file{std::string{path}};
  data.SerializeToOstream(&file);
}

void StoreTreeToProtobuf(MADAG dag, std::string_view path) {
  dag.AssertUA();
  Assert(dag.IsTree());

  Parsimony::data data;

  std::string newick;
  auto to_newick = [&](auto& self, MADAG::Node node) -> void {
    if (not node.IsLeaf()) {
      newick += '(';
    }
    size_t clade_idx = 0;
    for (auto clade : node.GetClades()) {
      Assert(clade.size() == 1);
      MADAG::Node i = (*clade.begin()).GetChild();
      if (i.IsLeaf()) {
        if (i.GetSampleId()) {
          newick += *i.GetSampleId();
        } else {
          newick += "unknown_leaf_";
          newick += std::to_string(i.GetId().value);
        }
      } else {
        self(self, i);
        newick += "inner_";
        newick += std::to_string(i.GetId().value);
      }
      if (++clade_idx < node.GetCladesCount()) {
        newick += ',';
      }
    }
    if (not node.IsLeaf()) {
      newick += ')';
    }
  };
  to_newick(to_newick, dag.GetRoot().GetFirstChild().GetChild());
  newick += ';';
  data.set_newick(newick);

  auto store_mutations = [](auto& self, MADAG::Edge edge, Parsimony::data& result,
                            std::string_view ref_seq) -> void {
    auto* proto = result.add_node_mutations();
    for (auto [pos, mut] : edge.GetEdgeMutations()) {
      auto* proto_mut = proto->add_mutation();
      proto_mut->set_position(static_cast<int32_t>(pos.value));
      proto_mut->set_ref_nuc(EncodeBase(ref_seq.at(pos.value - 1)));
      proto_mut->set_par_nuc(EncodeBase(mut.first));
      proto_mut->add_mut_nuc(EncodeBase(mut.second));
      proto_mut->set_chromosome("leaf_0");
    }
    for (MADAG::Edge child : edge.GetChild().GetChildren()) {
      self(self, child, result, ref_seq);
    }
  };
  store_mutations(store_mutations, dag.GetRoot().GetFirstChild(), data,
                  dag.GetReferenceSequence());

  std::ofstream file{std::string{path}};
  data.SerializeToOstream(&file);
}

static std::string EdgeMutationsToString(MADAG::Edge edge) {
  std::string result;
  size_t count = 0;
  for (auto [pos, muts] : edge.GetEdgeMutations()) {
    result += muts.first;
    result += std::to_string(pos.value);
    result += muts.second;
    result += ++count % 3 == 0 ? "\\n" : " ";
  }
  return result;
}

static std::string CompactGenomeToString(MADAG::Node node) {
  if (node.IsRoot()) {
    return "p";
  }
  std::string result = std::to_string(node.GetId().value);
  result += "\\n";
  size_t count = 0;
  for (auto [pos, base] : node.GetCompactGenome()) {
    result += std::to_string(pos.value);
    result += base;
    result += ++count % 3 == 0 ? "\\n" : " ";
  }
  return result;
}

void MADAGToDOT(MADAG dag, std::ostream& out) {
  out << "digraph {\n";
  out << "  forcelabels=true\n";
  out << "  nodesep=1.0\n";
  out << "  ranksep=2.0\n";
  out << "  ratio=1.0\n";
  for (MADAG::Edge edge : dag.GetEdges()) {
    auto [parent, child] = edge;
    out << "  \"" << CompactGenomeToString(parent) << "\" -> \""
        << CompactGenomeToString(child) << "\"";
    out << "[ xlabel=\"";
    out << EdgeMutationsToString(edge);
    out << "\" ]";
    out << "\n";
  }
  out << "}\n";
}

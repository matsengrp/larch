#include <fstream>
#include <vector>
#include <unordered_map>

#include <sys/stat.h>
#include <fcntl.h>

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/gzip_stream.h>

#include <range/v3/view/transform.hpp>

#include "nlohmann/json.hpp"

#include "dag_loader.hpp"
#include "dag.pb.h"
#include "parsimony.pb.h"
#include "newick.hpp"

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

MADAG LoadDAGFromProtobuf(std::string_view path) {
  ProtoDAG::data data;
  Parse(data, path);

  MADAG result;
  result.reference_sequence = data.reference_seq();

  for (auto& i : data.node_names()) {
    result.dag.AddNode({static_cast<size_t>(i.node_id())});
  }

  size_t edge_id = 0;
  for (auto& i : data.edges()) {
    result.dag.AddEdge({edge_id++}, {static_cast<size_t>(i.parent_node())},
                       {static_cast<size_t>(i.child_node())},
                       {static_cast<size_t>(i.parent_clade())});
  }

  result.dag.BuildConnections();

  result.edge_mutations.resize(result.dag.GetEdgesCount());
  edge_id = 0;
  for (auto& i : data.edges()) {
    EdgeMutations& muts = result.edge_mutations.at(edge_id++);
    for (auto& mut : i.edge_mutations()) {
      static const char decode[] = {'A', 'C', 'G', 'T'};
      Assert(mut.position() > 0);
      Assert(mut.mut_nuc().size() == 1);
      muts[{static_cast<size_t>(mut.position())}] = {decode[mut.par_nuc()],
                                                     decode[mut.mut_nuc().Get(0)]};
    }
  }

  result.dag.ReindexPreOrder();
  return result;
}

MADAG LoadTreeFromProtobuf(std::string_view path) {
  Parsimony::data data;
  Parse(data, path);

  MADAG result;

  size_t edge_id = 0;
  std::unordered_map<size_t, size_t> num_children;
  ParseNewick(
      data.newick(),
      [&result](size_t id, std::string label, std::optional<double> branch_length) {
        result.dag.AddNode({id});
        std::ignore = label;
        std::ignore = branch_length;
      },
      [&result, &edge_id, &num_children](size_t parent, size_t child) {
        result.dag.AddEdge({edge_id++}, {parent}, {child}, {num_children[parent]++});
      });
  result.dag.BuildConnections();

  result.edge_mutations.resize(result.dag.GetEdgesCount());

  size_t muts_idx = 0;
  for (Node node : result.dag.TraversePreOrder()) {
    const auto& pb_muts = data.node_mutations().Get(muts_idx++).mutation();
    if (node.IsRoot()) {
      continue;
    }
    auto& edge_muts = result.edge_mutations.at(node.GetSingleParent().GetId().value);
    for (auto i :
         pb_muts |
             ranges::views::transform(
                 [](auto& mut) -> std::pair<MutationPosition, std::pair<char, char>> {
                   static const char decode[] = {'A', 'C', 'G', 'T'};
                   Assert(mut.mut_nuc().size() == 1);
                   return {{static_cast<size_t>(mut.position())},
                           {decode[mut.par_nuc()], decode[mut.mut_nuc().Get(0)]}};
                 })) {
      edge_muts.insert(i);
    }
  }
  result.dag.ReindexPreOrder();
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
    // std::string par_nuc = mutation[1][0].get<std::string>();
    // Assert(par_nuc.size() == 1);
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

MADAG LoadDAGFromJson(std::string_view path) {
  nlohmann::json json = LoadJson(path);
  MADAG result;

  result.reference_sequence = json["refseq"][1];

  size_t id = 0;
  for ([[maybe_unused]] auto& i : json["nodes"]) {
    result.dag.AddNode({id++});
    size_t compact_genome_index = i[0];
    result.compact_genomes.push_back(GetCompactGenome(json, compact_genome_index));
  }
  id = 0;
  for (auto& i : json["edges"]) {
    result.dag.AddEdge({id++}, {i[0]}, {i[1]}, {i[2]});
  }
  result.dag.BuildConnections();
  std::map<NodeId, NodeId> index = result.dag.ReindexPreOrder();
  {
    std::vector<CompactGenome> reindexed_cgs;
    reindexed_cgs.reserve(result.compact_genomes.size());
    for (auto [prev, curr] : index) {
      reindexed_cgs.emplace_back(std::move(result.compact_genomes.at(curr.value)));
    }
    result.compact_genomes = std::move(reindexed_cgs);
  }
  return result;
}

void StoreDAGToProtobuf(const DAG& dag, std::string_view reference_sequence,
                        const std::vector<EdgeMutations>& edge_parent_mutations,
                        std::string_view path) {
  ProtoDAG::data data;

  data.set_reference_seq(std::string{reference_sequence});

  for (size_t i = 0; i < dag.GetNodesCount(); ++i) {
    auto* proto_node = data.add_node_names();
    proto_node->set_node_id(i);
  }

  for (Edge edge : dag.GetEdges()) {
    auto* proto_edge = data.add_edges();
    proto_edge->set_edge_id(edge.GetId().value);
    proto_edge->set_parent_node(edge.GetParentId().value);
    proto_edge->set_parent_clade(edge.GetChildId().value);
    proto_edge->set_child_node(edge.GetClade().value);
    for (auto [pos, nucs] : edge_parent_mutations.at(edge.GetId().value)) {
      auto* mut = proto_edge->add_edge_mutations();
      mut->set_position(pos.value);
      switch (nucs.first) {
        case 'A':
          mut->set_par_nuc(0);
          break;
        case 'C':
          mut->set_par_nuc(1);
          break;
        case 'G':
          mut->set_par_nuc(2);
          break;
        case 'T':
          mut->set_par_nuc(3);
          break;
        default:
          Fail("Invalid base");
      };
      switch (nucs.second) {
        case 'A':
          mut->add_mut_nuc(0);
          break;
        case 'C':
          mut->add_mut_nuc(1);
          break;
        case 'G':
          mut->add_mut_nuc(2);
          break;
        case 'T':
          mut->add_mut_nuc(3);
          break;
        default:
          Fail("Invalid base");
      };
    }
  }

  std::ofstream file{std::string{path}};
  data.SerializeToOstream(&file);
}

static std::string EdgeMutationsToString(Edge edge, const MADAG& dag) {
  std::string result;
  size_t count = 0;
  for (auto [pos, muts] : dag.edge_mutations.at(edge.GetId().value)) {
    result += muts.first;
    result += std::to_string(pos.value);
    result += muts.second;
    result += ++count % 3 == 0 ? "\\n" : " ";
  }
  return result;
}

static std::string CompactGenomeToString(
    Node node, const std::vector<CompactGenome>& compact_genomes) {
  if (node.IsRoot()) {
    return "p";
  }
  std::string result = std::to_string(node.GetId().value);
  result += "\\n";
  size_t count = 0;
  for (auto [pos, base] : compact_genomes.at(node.GetId().value)) {
    result += std::to_string(pos.value);
    result += base;
    result += ++count % 3 == 0 ? "\\n" : " ";
  }
  return result;
}

void MADAGToDOT(const MADAG& dag, std::ostream& out) {
  out << "digraph {\n";
  out << "  forcelabels=true\n";
  out << "  nodesep=1.0\n";
  out << "  ranksep=2.0\n";
  out << "  ratio=1.0\n";
  for (Edge edge : dag.dag.GetEdges()) {
    auto [parent, child] = edge;
    if (dag.compact_genomes.empty()) {
      out << parent.GetId().value << " -> " << child.GetId().value;
    } else {
      out << "  \"" << CompactGenomeToString(parent, dag.compact_genomes) << "\" -> \""
          << CompactGenomeToString(child, dag.compact_genomes) << "\"";
    }
    if (not dag.edge_mutations.empty()) {
      out << "[ xlabel=\"";
      out << EdgeMutationsToString(edge, dag);
      out << "\" ]";
    }
    out << "\n";
  }
  out << "}\n";
}

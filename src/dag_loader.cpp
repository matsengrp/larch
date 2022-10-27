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

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>

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

  MADAG result{data.reference_seq()};

  for (auto& i : data.node_names()) {
    result.AddNode({static_cast<size_t>(i.node_id())});
  }

  size_t edge_id = 0;
  for (auto& i : data.edges()) {
    result.AddEdge({edge_id++}, {static_cast<size_t>(i.parent_node())},
                   {static_cast<size_t>(i.child_node())},
                   {static_cast<size_t>(i.parent_clade())});
  }

  result.BuildConnections();

  std::vector<EdgeMutations> edge_mutations;
  edge_mutations.resize(result.GetDAG().GetEdgesCount());
  edge_id = 0;
  for (auto& i : data.edges()) {
    EdgeMutations& muts = edge_mutations.at(edge_id++);
    for (auto& mut : i.edge_mutations()) {
      static const char decode[] = {'A', 'C', 'G', 'T'};
      Assert(mut.position() > 0);
      Assert(mut.mut_nuc().size() == 1);
      muts[{static_cast<size_t>(mut.position())}] = {decode[mut.par_nuc()],
                                                     decode[mut.mut_nuc().Get(0)]};
    }
  }
  result.SetEdgeMutations(std::move(edge_mutations));
  result.AssertUA();
  return result;
}

MADAG LoadTreeFromProtobuf(std::string_view path, std::string_view reference_sequence) {
  Parsimony::data data;
  Parse(data, path);

  MADAG result{reference_sequence};

  std::unordered_map<size_t, size_t> num_children;
  std::map<size_t, std::optional<std::string>> seq_ids;
  ParseNewick(
      data.newick(),
      [&result, &seq_ids](size_t node_id, std::string label, std::optional<double>) {
        seq_ids[node_id] = label;
      },
      [&result, &num_children](size_t parent, size_t child) {
        result.AddEdge({child}, {parent}, {child}, {num_children[parent]++});
      });
  result.InitializeNodes(result.GetDAG().GetEdgesCount() + 1);
  result.BuildConnections();

  for (auto node : result.GetDAG().GetNodes()) {
    if (node.IsLeaf()) {
      node.SetSampleId(std::move(seq_ids[node.GetId().value]));
    }
  }

  std::vector<EdgeMutations> edge_mutations;
  edge_mutations.resize(result.GetDAG().GetEdgesCount());
  Assert(static_cast<size_t>(data.node_mutations_size()) ==
         result.GetDAG().GetNodesCount());

  size_t muts_idx = 0;
  for (auto iter : result.GetDAG().TraversePreOrder()) {
    const auto& pb_muts = data.node_mutations().Get(muts_idx++).mutation();
    auto& edge_muts = edge_mutations.at(iter.GetEdge().GetId().value);
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
  result.SetEdgeMutations(std::move(edge_mutations));

  // add UA node to tree, with mutations corresponding to original tree root node mutations (if any)
  const auto& pb_muts = data.node_mutations().Get(0).mutation();
  EdgeMutations edge_muts;
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

  result.AddUA(std::move(edge_muts));

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

MADAG LoadDAGFromJson(std::string_view path) {
  nlohmann::json json = LoadJson(path);
  MADAG result{std::string{json["refseq"][1]}};

  size_t id = 0;
  for ([[maybe_unused]] auto& i : json["nodes"]) {
    result.AddNode({id++});
    size_t compact_genome_index = i[0];
    result.AppendCompactGenome(GetCompactGenome(json, compact_genome_index));
  }
  id = 0;
  for (auto& i : json["edges"]) {
    result.AddEdge({id++}, {i[0]}, {i[1]}, {i[2]});
  }
  result.BuildConnections();
  result.AssertUA();
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
    proto_edge->set_child_node(edge.GetChildId().value);
    proto_edge->set_parent_clade(edge.GetClade().value);
    for (auto [pos, nucs] : edge_parent_mutations.at(edge.GetId().value)) {
      auto* proto_mut = proto_edge->add_edge_mutations();
      proto_mut->set_position(pos.value);
      proto_mut->set_par_nuc(EncodeBase(nucs.first));
      proto_mut->add_mut_nuc(EncodeBase(nucs.second));
    }
  }

  std::ofstream file{std::string{path}};
  data.SerializeToOstream(&file);
}

void StoreDAGToProtobuf(const MADAG& dag, std::string_view path) {
  Assert(not dag.GetEdgeMutations().empty());
  dag.AssertUA();
  StoreDAGToProtobuf(dag.GetDAG(), dag.GetReferenceSequence(), dag.GetEdgeMutations(),
                     path);
}

void StoreTreeToProtobuf(const DAG& dag, std::string_view reference_sequence,
                         const std::vector<EdgeMutations>& edge_parent_mutations,
                         std::string_view path) {
  Assert(dag.IsTree());

  Parsimony::data data;

  std::string newick;
  auto to_newick = [&](auto& self, Node node) -> void {
    if (not node.IsLeaf()) {
      newick += '(';
    }
    size_t clade_idx = 0;
    for (auto clade : node.GetClades()) {
      Assert(clade.size() == 1);
      Node i = (*clade.begin()).GetChild();
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

  for (Edge edge : dag.TraversePreOrder(dag.GetRoot().GetFirstChild().GetChild())) {
    Assert(not edge.IsRoot());

    auto* proto = data.add_node_mutations();
    for (auto [pos, mut] : edge_parent_mutations.at(edge.GetId().value)) {
      auto* proto_mut = proto->add_mutation();
      proto_mut->set_position(static_cast<int32_t>(pos.value));
      proto_mut->set_ref_nuc(EncodeBase(reference_sequence.at(pos.value - 1)));
      proto_mut->set_par_nuc(EncodeBase(mut.first));
      proto_mut->add_mut_nuc(EncodeBase(mut.second));
      proto_mut->set_chromosome("leaf_0");
    }
  }

  std::ofstream file{std::string{path}};
  data.SerializeToOstream(&file);
}

void StoreTreeToProtobuf(const MADAG& dag, std::string_view path) {
  dag.AssertUA();
  Assert(not dag.GetEdgeMutations().empty());
  StoreTreeToProtobuf(dag.GetDAG(), dag.GetReferenceSequence(), dag.GetEdgeMutations(),
                      path);
}

static std::string EdgeMutationsToString(Edge edge, const MADAG& dag) {
  std::string result;
  size_t count = 0;
  for (auto [pos, muts] : dag.GetEdgeMutations(edge)) {
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
  for (Edge edge : dag.GetDAG().GetEdges()) {
    auto [parent, child] = edge;
    if (dag.GetCompactGenomes().empty()) {
      out << parent.GetId().value << " -> " << child.GetId().value;
    } else {
      out << "  \"" << CompactGenomeToString(parent, dag.GetCompactGenomes())
          << "\" -> \"" << CompactGenomeToString(child, dag.GetCompactGenomes())
          << "\"";
    }
    if (not dag.GetEdgeMutations().empty()) {
      out << "[ xlabel=\"";
      out << EdgeMutationsToString(edge, dag);
      out << "\" ]";
    }
    out << "\n";
  }
  out << "}\n";
}

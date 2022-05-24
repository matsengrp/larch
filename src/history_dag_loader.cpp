#include <fstream>
#include <vector>
#include <unordered_map>

#include <range/v3/view/transform.hpp>

#include "nlohmann/json.hpp"

#include "history_dag_loader.hpp"
#include "zlib_stream.hpp"
#include "dag.pb.h"
#include "parsimony.pb.h"
#include "newick.hpp"

HistoryDAG LoadHistoryDAGFromProtobufGZ(std::string_view path, std::string& ref_seq,
                                        std::vector<Mutations>& mutations) {
  std::ifstream in_compressed{std::string{path}};
  assert(in_compressed);
  zlib::ZStringBuf zbuf{in_compressed, 1, 128 * 1024 * 1024};
  std::istream in{&zbuf};
  DAG::data data;
  data.ParseFromIstream(&in);

  ref_seq = data.reference_seq();
  HistoryDAG dag;

  for (auto& i : data.node_names()) {
    dag.AddNode({static_cast<size_t>(i.node_id())});
  }

  size_t edge_id = 0;
  for (auto& i : data.edges()) {
    dag.AddEdge({edge_id++}, {static_cast<size_t>(i.parent_node())},
                {static_cast<size_t>(i.child_node())},
                {static_cast<size_t>(i.parent_clade())});
  }

  dag.BuildConnections();

  mutations.resize(dag.GetEdges().size());
  edge_id = 0;
  for (auto& i : data.edges()) {
    Mutations& cg = mutations.at(edge_id++);
    for (auto& mut : i.edge_mutations()) {
      static const char decode[] = {'A', 'C', 'G', 'T'};
      assert(mut.mut_nuc().size() == 1);
      cg[{static_cast<size_t>(mut.position())}] = decode[mut.mut_nuc().Get(0)];
    }
  }

  return dag;
}

static HistoryDAG LoadTreeFromProtobuf(std::istream& in,
                                       std::vector<Mutations>& mutations) {
  Parsimony::data data;
  data.ParseFromIstream(&in);

  HistoryDAG dag;

  size_t edge_id = 0;
  std::unordered_map<size_t, size_t> num_children;
  ParseNewick(
      data.newick(),
      [&dag](size_t id, std::string label, std::optional<double> branch_length) {
        dag.AddNode({id});
        std::ignore = label;
        std::ignore = branch_length;
      },
      [&dag, &edge_id, &num_children](size_t parent, size_t child) {
        dag.AddEdge({edge_id++}, {parent}, {child}, {num_children[parent]++});
      });
  dag.BuildConnections();

  mutations.resize(dag.GetEdges().size());

  size_t muts_idx = 0;
  for (MutableNode node : dag.TraversePreOrder()) {
    const auto& pb_muts = data.node_mutations().Get(muts_idx++).mutation();
    if (node.IsRoot()) {
      continue;
    }
    auto& edge_muts = mutations.at(node.GetSingleParent().GetId().value);
    for (auto i :
         pb_muts | ranges::view::transform([](auto& mut)
                                               -> std::pair<MutationPosition, char> {
           static const char decode[] = {'A', 'C', 'G', 'T'};
           assert(mut.mut_nuc().size() == 1);
           return {{static_cast<size_t>(mut.position())}, decode[mut.mut_nuc().Get(0)]};
         })) {
      edge_muts.insert(i);
    }
  }

  return dag;
}

static bool IsGzipped(std::string_view path) {
  std::ifstream in{std::string{path}};
  assert(in);
  unsigned char header[2];
  in >> header[0] >> header[1];
  return header[0] == 0x1f and header[1] == 0x8b;
}

HistoryDAG LoadTreeFromProtobufGZ(std::string_view path,
                                  std::vector<Mutations>& mutations) {
  if (IsGzipped(path)) {
    std::ifstream in_compressed{std::string{path}};
    assert(in_compressed);
    zlib::ZStringBuf zbuf{in_compressed, 1, 128 * 1024 * 1024};
    std::istream in{&zbuf};
    return LoadTreeFromProtobuf(in, mutations);
  } else {
    std::ifstream in{std::string{path}};
    return LoadTreeFromProtobuf(in, mutations);
  }
}

[[nodiscard]] std::string LoadRefseqFromJsonGZ(std::string_view path) {
  std::ifstream in_compressed{std::string{path}};
  assert(in_compressed);
  zlib::ZStringBuf zbuf{in_compressed, 1, 128 * 1024 * 1024};
  std::istream in{&zbuf};
  nlohmann::json json;
  in >> json;
  return json["refseq"][1];
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

HistoryDAG LoadHistoryDAGFromJsonGZ(std::string_view path, std::string& refseq) {
  std::ifstream in_compressed{std::string{path}};
  assert(in_compressed);
  zlib::ZStringBuf zbuf{in_compressed, 1, 128 * 1024 * 1024};
  std::istream in{&zbuf};
  nlohmann::json json;
  in >> json;
  HistoryDAG result;

  refseq = json["refseq"][1];

  size_t id = 0;
  for ([[maybe_unused]] auto& i : json["nodes"]) {
    result.AddNode({id++});
  }
  id = 0;
  for (auto& i : json["edges"]) {
    result.AddEdge({id++}, {i[0]}, {i[1]}, {i[2]});
  }
  result.BuildConnections();
  return result;
}

// static CompactGenome GetCompactGenome(const nlohmann::json& json,
//                                       size_t compact_genome_index) {
//   CompactGenome result;
//   for (auto& mutation : json["compact_genomes"][compact_genome_index]) {
//     size_t position = mutation[0];
//     std::string new_base = mutation[1][1].get<std::string>();
//     result[position] = new_base.at(0);
//   }
//   return result;
// }

// static LeafSet GetLeafSet(const nlohmann::json& json, size_t node_index) {
//   LeafSet result;
//   for (auto& clade : json["nodes"][node_index][1]) {
//     std::set<CompactGenome*> leafs;
//     for (size_t leaf : clade) {
//       leafs.insert(GetCompactGenome(json, leaf));
//     }
//     result.insert(leafs);
//   }
//   return result;
// }

// std::vector<NodeLabel> LoadLabelsJsonGZ(std::string_view path) {
//   std::ifstream in_compressed{std::string{path}};
//   assert(in_compressed);
//   zlib::ZStringBuf zbuf{in_compressed, 1, 128 * 1024 * 1024};
//   std::istream in{&zbuf};
//   nlohmann::json json;
//   in >> json;

//   std::vector<NodeLabel> result;
//   size_t node_index = 0;
//   for (auto& node : json["nodes"]) {
//     size_t compact_genome_index = node[0];
//     result.push_back(
//         {GetCompactGenome(json, compact_genome_index), GetLeafSet(json,
//         node_index++)});
//   }
//   return result;
// }

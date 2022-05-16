#include <fstream>
#include <vector>
#include <unordered_map>

#include "nlohmann/json.hpp"

#include "history_dag_loader.hpp"
#include "zlib_stream.hpp"
#include "dag.pb.h"
#include "newick.hpp"

HistoryDAG LoadHistoryDAGFromProtobufGZ(const std::string& path, std::vector<CompactGenome>& mutations) {
    std::ifstream in_compressed(path.c_str());
    assert(in_compressed);
    zlib::ZStringBuf zbuf(in_compressed, 1, 128 * 1024 * 1024);
    std::istream in(&zbuf);
    DAG::data data;
    data.ParseFromIstream(&in);
    HistoryDAG dag;

    for (auto& i : data.node_names()) {
        dag.AddNode({static_cast<size_t>(i.node_id())});
    }

    size_t edge_id = 0;
    for (auto& i : data.edges()) {
        dag.AddEdge(
            {edge_id++},
            {static_cast<size_t>(i.parent_node())},
            {static_cast<size_t>(i.child_node())},
            {static_cast<size_t>(i.parent_clade())});
    }

    dag.BuildConnections();

    mutations.resize(dag.GetEdges().size());
    edge_id = 0;
    for (auto& i : data.edges()) {
        CompactGenome& cg = mutations.at(edge_id++);
        for (auto& mut : i.edge_mutations()) {
            static const char decode[] = {'A', 'C', 'G', 'T'};
            assert(mut.mut_nuc().size() == 1);
            cg[mut.position()] = decode[mut.mut_nuc().Get(0)];
        }
    }

    return dag;
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

HistoryDAG LoadHistoryDAGFromJsonGZ(const std::string& path, std::string& refseq) {
        
    std::ifstream in_compressed(path);
    assert(in_compressed);
    zlib::ZStringBuf zbuf(in_compressed, 1, 128 * 1024 * 1024);
    std::istream in(&zbuf);
    nlohmann::json json;
    in >> json;
    HistoryDAG result;

    refseq = json["refseq"][1];

    size_t id = 0;
    for ([[maybe_unused]] auto& i : json["nodes"]) result.AddNode({id++});
    id = 0;
    for (auto& i : json["edges"]) {
        MutableEdge edge = result.AddEdge({id++}, {i[0]}, {i[1]}, {i[2]});
        std::ignore = edge;
        // size_t child_id = i[1];
        // size_t cg_idx = json["nodes"][child_id][0];
        // std::map<MutationPosition, char>& cg = GetOrInsert(genomes, child_id);
        // for (auto& j : json["compact_genomes"][cg_idx]) {
        //     cg.push_back({
        //         {static_cast<size_t>(j[0])},
        //         j[1][0].get<std::string>()[0],
        //         j[1][1].get<std::string>()[0]
        //     });
        // }
        // edge.SetMutations(ranges::views::all(cg));
    }
    result.BuildConnections();

    // for (MutableEdge i : result.TraversePreOrder()) {
    //     if (i.IsRoot()) continue;
    //     auto parent_muts = i.GetParent().GetFirstParent().GetMutations();
    //     std::set<Mutation> mut_set;
    //     mut_set.insert(parent_muts.begin(), parent_muts.end());
    //     std::vector<Mutation> muts;
    //     for (auto j : i.GetMutations()) {
    //         if (mut_set.find(j) == mut_set.end()) muts.push_back(j);
    //     }
    //     i.SetMutations(muts | ranges::views::all);
    // }

    return result;
}

static CompactGenome GetCompactGenome(const nlohmann::json& json, size_t compact_genome_index) {
	CompactGenome result;
	for (auto& mutation : json["compact_genomes"][compact_genome_index]) {
		size_t position = mutation[0];
		std::string new_base = mutation[1][1].get<std::string>();
		result[position] = new_base.at(0);
	}
	return result;
}

static LeafSet GetLeafSet(const nlohmann::json& json, size_t node_index) {
	LeafSet result;
	for (auto& clade : json["nodes"][node_index][1]) {
		std::set<CompactGenome> leafs;
		for (size_t leaf : clade) {
			leafs.insert(GetCompactGenome(json, leaf));
		}
		result.insert(leafs);
	}
	return result;
}

std::vector<NodeLabel> LoadLabelsJsonGZ(const std::string& path) {
    std::ifstream in_compressed(path);
    assert(in_compressed);
    zlib::ZStringBuf zbuf(in_compressed, 1, 128 * 1024 * 1024);
    std::istream in(&zbuf);
    nlohmann::json json;
    in >> json;

    std::vector<NodeLabel> result;
    size_t node_index = 0;
	for (auto& node : json["nodes"]) {
		size_t compact_genome_index = node[0];
		result.push_back({GetCompactGenome(json, compact_genome_index),
			GetLeafSet(json, node_index++)});
	}
    return result;
}

#include "subtree_weight.hpp"
#include "parsimony_score.hpp"
#include "merge.hpp"

#include <iostream>
#include <fstream>
#include <string_view>
#include <unordered_set>
#include <algorithm>

#include "test_common.hpp"

#include "dag_loader.hpp"

bool compare_treedags(MADAG& dag1, MADAG& dag2) {
    if (dag1.GetReferenceSequence() != dag2.GetReferenceSequence()) {
        return false;
    }
    if (dag1.GetDAG().GetNodesCount() != dag2.GetDAG().GetNodesCount()) {
        return false;
    }
    if (dag1.GetDAG().GetEdgesCount() != dag2.GetDAG().GetEdgesCount()) {
        return false;
    }

    if (not dag1.GetReferenceSequence().empty()) {
        if (dag1.GetCompactGenomes().empty()) {
            dag1.GetCompactGenomes() = dag1.ComputeCompactGenomes(dag1.GetReferenceSequence());
            dag2.GetCompactGenomes() = dag2.ComputeCompactGenomes(dag2.GetReferenceSequence());
        }
        std::unordered_set<CompactGenome> dag1_cgs, dag2_cgs;
        for (auto &cg : dag1.GetCompactGenomes()) {
             dag1_cgs.emplace(cg.Copy());
        }
        for (auto &cg : dag2.GetCompactGenomes()) {
             dag2_cgs.emplace(cg.Copy());
        }
        if (dag1_cgs != dag2_cgs) {
            return false;
        }
    }

    if (not (dag1.GetEdgeMutations().empty() or dag2.GetEdgeMutations().empty())) {
        std::vector<EdgeMutations> &dag1_ems = dag1.GetEdgeMutations();
        std::vector<EdgeMutations> &dag2_ems = dag2.GetEdgeMutations();
        for (auto &em : dag1_ems) {
            if (std::count(dag2_ems.begin(), dag2_ems.end(), em) != std::count(dag1_ems.begin(), dag1_ems.end(), em)) {
                return false;
            }
        }
    } else if (not (dag1.GetEdgeMutations().empty() and dag2.GetEdgeMutations().empty())) {
        return false;
    }
    return true;
}

/*
std::unordered_map map_between_leaf_nodes(MADAG dag1, MADAG dag2) {
    std::vector<CompactGenome> &dag1_cgs = dag1.GetCompactGenomes();
    std::vector<CompactGenome> &dag2_cgs = dag2.GetCompactGenomes();

    if (dag1_cgs.size() > 0 and dag2_cgs.size() > 0) {
        std::unordered_map<NodeId, NodeId> map_between_leaf_nodes;

        for (auto node1 : dag1.GetDAG().GetNodes()) {
            if (node1.IsLeaf()) {
                CompactGenome cg1 = dag1_cgs.at(node1.GetId().value).Copy();
                for (auto node2 : dag2.GetDAG().GetNodes()) {
                     if (node2.IsLeaf() and 
                         (cg1 == dag2_cgs.at(node2.GetId().value).Copy())) {
                         map_between_leaf_nodes[node1.GetId()] = node2.GetId();
                     }
                }
            }
        }
    }
    return map_between_leaf_nodes;
}
*/
void assign_sample_ids(MADAG &dag1, MADAG &dag2) {
    std::vector<CompactGenome> &dag1_cgs = dag1.GetCompactGenomes();
    std::vector<CompactGenome> &dag2_cgs = dag2.GetCompactGenomes();

    if (dag1_cgs.size() > 0 and dag2_cgs.size() > 0) {
        for (auto node1 : dag1.GetDAG().GetNodes()) {
            if (node1.IsLeaf()) {
                CompactGenome cg1 = dag1_cgs.at(node1.GetId().value).Copy();
                for (auto node2 : dag2.GetDAG().GetNodes()) {
                     if (node2.IsLeaf()  and 
                         (cg1 == dag2_cgs.at(node2.GetId().value).Copy())) {
                         node2.SetSampleId(node1.GetSampleId());
                     }
                }
            }
        }
    }
}

static void test_write_protobuf() {
    std::string_view path = "data/check_parsimony_protobuf/example_tree.pb";
    MADAG treedag = LoadTreeFromProtobuf(path);
    std::fstream file;
    std::string refseq, filename;
    filename = "data/check_parsimony_protobuf/refseq.fasta";
    file.open(filename);
    while (file >> refseq){}
    treedag.GetReferenceSequence() = refseq;
    std::cout << "loaded...\n" << std::flush;

    SubtreeWeight<ParsimonyScore> weight{treedag};
    MADAG sample_tree = weight.SampleTree({});
    std::cout << "sampled the dag...\n" << std::flush;

    std::cout << "comparing tree to its sample\n" << std::flush;
    //StoreTreeToProtobuf(sample_tree, "/home/wdumm/larch/test_write_protobuf.pb");
    compare_treedags(treedag, sample_tree);

    sample_tree.GetCompactGenomes() = sample_tree.ComputeCompactGenomes(sample_tree.GetReferenceSequence());
    treedag.GetCompactGenomes() = treedag.ComputeCompactGenomes(treedag.GetReferenceSequence());
    std::cout << "comparing tree to its sample, with computed compact genomes\n" << std::flush;
    compare_treedags(treedag, sample_tree);

    assign_sample_ids(treedag, sample_tree);
    for (auto node : sample_tree.GetDAG().GetNodes()) {
        if (node.IsLeaf()) {
            std::cout << node.GetSampleId().value() << ", ";
        }
    }

    //StoreTreeToProtobuf(sample_tree, "/home/wdumm/larch/test_write_protobuf.pb");

    treedag.GetEdgeMutations() = treedag.ComputeEdgeMutations(treedag.GetReferenceSequence());
    std::cout << "comparing recomputed edge mutations to original\n" << std::flush;
    compare_treedags(treedag, sample_tree);

    Merge merge{treedag.GetReferenceSequence()};
    merge.AddDAGs({treedag, sample_tree});
    merge.GetResult().GetEdgeMutations() = merge.ComputeResultEdgeMutations();

    std::cout << "comparing original treedag to trivial merge\n" << std::flush;
    compare_treedags(treedag, merge.GetResult());

    for (auto node : treedag.GetDAG().GetNodes()) {
        if (node.IsLeaf()) {
            std::cout << node.GetSampleId().value() << ", ";
        }
    }
    std::cout << "\n";
    StoreDAGToProtobuf(treedag.GetDAG(),
                    treedag.GetReferenceSequence(),
                    treedag.GetEdgeMutations(),
                    "data/check_parsimony_protobuf/example_tree_in_dag_format.pb");

    std::cout << "wrote dag to protobuf...\n" << std::flush;
}

[[maybe_unused]] static const auto test_added_write =
    add_test({test_write_protobuf, "test write protobuf"});

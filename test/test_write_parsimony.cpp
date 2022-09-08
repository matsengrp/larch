#include "subtree_weight.hpp"
#include "parsimony_score.hpp"
#include "merge.hpp"

#include <iostream>
#include <fstream>
#include <string_view>

#include "test_common.hpp"

#include "dag_loader.hpp"

static void test_write_protobuf() {
    std::string_view path = "data/check_parsimony_protobuf/example_tree.pb";
    MADAG treedag = LoadTreeFromProtobuf(path);
    std::fstream file;
    std::string refseq, filename;
    filename = "data/check_parsimony_protobuf/refseq.fasta";
    file.open(filename);
    while (file >> refseq){}
    treedag.GetReferenceSequence() = refseq;
    std::cout << "loaded..." << std::flush;

    size_t counter = 0;
    for (auto node : treedag.GetDAG().GetNodes()) {
        if (node.IsLeaf()) {
            counter++;
        }
    }

    std::cout << "tree has " << counter << " leaf nodes" << std::flush;
    treedag.GetCompactGenomes() = treedag.ComputeCompactGenomes(treedag.GetReferenceSequence());
    treedag.GetEdgeMutations() = treedag.ComputeEdgeMutations(treedag.GetReferenceSequence());
    
    /* std::cout << "computed compact genomes on original treedag...\n" << std::flush; */

    /* SubtreeWeight<ParsimonyScore> weight{treedag}; */
    /* MADAG sample_tree = weight.SampleTree({}); */
    /* std::cout << "sampled...\n" << std::flush; */
    /* auto cgs1_ = sample_tree.ComputeCompactGenomes(sample_tree.GetReferenceSequence()); */
    /* std::cout << "computed compact genomes on sampled tree...\n" << std::flush; */

    /* Merge merge{treedag.GetReferenceSequence()}; */
    /* merge.AddDAGs({treedag, sample_tree}); */


    StoreTreeToProtobuf(treedag, "/home/wdumm/larch/test_write_protobuf.pb");
}

[[maybe_unused]] static const auto test_added_write =
    add_test({test_write_protobuf, "test write protobuf"});

#include "subtree_weight.hpp"
#include "parsimony_score.hpp"
#include "merge.hpp"

#include <iostream>
#include <fstream>
#include <string_view>

#include "test_common.hpp"

#include "dag_loader.hpp"

bool compare_treedags(MADAG& dag1, MADAG& dag2) {
    assert (dag1.GetReferenceSequence() == dag2.GetReferenceSequence());
    assert (dag1.GetDAG().GetNodesCount() == dag2.GetDAG().GetNodesCount());
    assert (dag1.GetDAG().GetEdgesCount() == dag2.GetDAG().GetEdgesCount());
    /* assert (not (not dag1.GetCompactGenomes().empty() and dag2.GetCompactGenomes().empty())); */
    /* assert (not (not dag1.GetEdgeMutations().empty() and dag2.GetEdgeMutations().empty())); */
    assert(dag1.GetCompactGenomes() == dag2.GetCompactGenomes());
    assert(dag1.GetEdgeMutations() == dag2.GetEdgeMutations());
    return true;
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
    std::cout << "loaded..." << std::flush;

    SubtreeWeight<ParsimonyScore> weight{treedag};
    MADAG sample_tree = weight.SampleTree({});

    std::cout << "comparing tree to its sample\n" << std::flush;
    StoreTreeToProtobuf(sample_tree, "/home/wdumm/larch/test_write_protobuf.pb");
    compare_treedags(treedag, sample_tree);
    sample_tree.GetCompactGenomes() = sample_tree.ComputeCompactGenomes(sample_tree.GetReferenceSequence());
    treedag.GetCompactGenomes() = treedag.ComputeCompactGenomes(treedag.GetReferenceSequence());
    std::cout << "comparing tree to its sample, with computed compact genomes\n" << std::flush;
    compare_treedags(treedag, sample_tree);

    treedag.GetEdgeMutations() = treedag.ComputeEdgeMutations(treedag.GetReferenceSequence());
    std::cout << "comparing recomputed edge mutations to original\n" << std::flush;
    compare_treedags(treedag, sample_tree);

    /* Merge merge{treedag.GetReferenceSequence()}; */
    /* merge.AddDAGs({treedag, sample_tree}); */
    /* merge.GetResult().GetEdgeMutations() = merge.ComputeResultEdgeMutations(); */


}

[[maybe_unused]] static const auto test_added_write =
    add_test({test_write_protobuf, "test write protobuf"});

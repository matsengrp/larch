#include "subtree_weight.hpp"
#include "parsimony_score.hpp"

#include <iostream>
#include <fstream>
#include <string_view>

#include "test_common.hpp"

#include "dag_loader.hpp"

static void test_write_protobuf() {
    std::string_view path = "data/check_parsimony_protobuf/example_tree.pb";
    MADAG treedag = LoadTreeFromProtobuf(path);
    std::cout << "loaded..." << std::flush;

    std::fstream file;
    std::string refseq, filename;
    filename = "data/check_parsimony_protobuf/refseq.fasta";
    file.open(filename);
    while (file >> refseq){}

    treedag.GetReferenceSequence() = refseq;
    std::cout << refseq << std::flush;
    StoreTreeToProtobuf(treedag, "/home/wdumm/larch/test_write_protobuf.pb");
}

[[maybe_unused]] static const auto test_added_write =
    add_test({test_write_protobuf, "test write protobuf"});

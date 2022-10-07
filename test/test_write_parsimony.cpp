#include "subtree_weight.hpp"
#include "parsimony_score_binary.hpp"
#include "merge.hpp"

#include <iostream>
#include <fstream>
#include <string_view>
#include <unordered_set>
#include <algorithm>

#include "test_common.hpp"

#include "dag_loader.hpp"

bool compare_treedags(MADAG &dag1, MADAG &dag2) {
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
      dag1.RecomputeCompactGenomes();
      dag2.RecomputeCompactGenomes();
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

  if (not(dag1.GetEdgeMutations().empty() or dag2.GetEdgeMutations().empty())) {
    const std::vector<EdgeMutations> &dag1_ems = dag1.GetEdgeMutations();
    const std::vector<EdgeMutations> &dag2_ems = dag2.GetEdgeMutations();

    for (auto &em : dag1_ems) {
      if (std::count(dag2_ems.begin(), dag2_ems.end(), em) !=
          std::count(dag1_ems.begin(), dag1_ems.end(), em)) {
        return false;
      }
    }
  } else if (not(dag1.GetEdgeMutations().empty() and dag2.GetEdgeMutations().empty())) {
    return false;
  }
  return true;
}

static void test_write_protobuf() {
  std::string_view path = "data/check_parsimony_protobuf/example_tree.pb";
  std::fstream file;
  std::string refseq, filename;
  filename = "data/check_parsimony_protobuf/refseq.fasta";
  file.open(filename);
  while (file >> refseq) {
  }
  MADAG treedag = LoadTreeFromProtobuf(path, refseq);

  SubtreeWeight<ParsimonyScore> weight{treedag};
  MADAG sample_tree = weight.SampleTree({}).first;

  StoreTreeToProtobuf(sample_tree, "test_write_protobuf.pb");
  compare_treedags(treedag, sample_tree);

  sample_tree.RecomputeCompactGenomes();
  treedag.RecomputeCompactGenomes();
  compare_treedags(treedag, sample_tree);

  treedag.RecomputeEdgeMutations();
  compare_treedags(treedag, sample_tree);

  Merge merge{treedag.GetReferenceSequence()};
  merge.AddDAGs({treedag, sample_tree});
  merge.ComputeResultEdgeMutations();

  compare_treedags(treedag, merge.GetResult());
}

[[maybe_unused]] static const auto test_added_write =
    add_test({test_write_protobuf, "test write protobuf"});

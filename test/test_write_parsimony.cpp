#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score_binary.hpp"
#include "larch/merge/merge.hpp"

#include <iostream>
#include <fstream>
#include <string_view>
#include <unordered_set>
#include <algorithm>

#include "test_common.hpp"

#include "larch/dag_loader.hpp"

template <typename DAG1, typename DAG2>
bool compare_treedags(DAG1 dag1, DAG2 dag2) {
  if (dag1.GetReferenceSequence() != dag2.GetReferenceSequence()) {
    return false;
  }
  if (dag1.GetNodesCount() != dag2.GetNodesCount()) {
    return false;
  }
  if (dag1.GetEdgesCount() != dag2.GetEdgesCount()) {
    return false;
  }

  std::unordered_set<CompactGenome> dag1_cgs, dag2_cgs;
  for (auto node : dag1.GetNodes()) {
    dag1_cgs.emplace(node.GetCompactGenome().Copy());
  }
  for (auto node : dag2.GetNodes()) {
    dag2_cgs.emplace(node.GetCompactGenome().Copy());
  }
  if (dag1_cgs != dag2_cgs) {
    return false;
  }

  std::vector<EdgeMutations> dag1_ems;
  std::vector<EdgeMutations> dag2_ems;

  for (auto edge : dag1.GetEdges()) {
    dag1_ems.emplace_back(edge.GetEdgeMutations().Copy());
  }
  for (auto edge : dag2.GetEdges()) {
    dag2_ems.emplace_back(edge.GetEdgeMutations().Copy());
  }

  if (not(dag1_ems.empty() or dag2_ems.empty())) {
    for (auto &em : dag1_ems) {
      if (std::count(dag2_ems.begin(), dag2_ems.end(), em) !=
          std::count(dag1_ems.begin(), dag1_ems.end(), em)) {
        return false;
      }
    }
  } else if (not(dag1_ems.empty() and dag2_ems.empty())) {
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
  MADAGStorage treedag = LoadTreeFromProtobuf(path, refseq);

  SubtreeWeight<BinaryParsimonyScore, MADAG> weight{treedag.View()};
  auto sample_tree = weight.SampleTree({});

  StoreTreeToProtobuf(sample_tree.View(), "test_write_protobuf.pb");
  compare_treedags(treedag.View(), sample_tree.View());

  sample_tree.View().RecomputeCompactGenomes();
  treedag.View().RecomputeCompactGenomes();
  compare_treedags(treedag.View(), sample_tree.View());

  treedag.View().RecomputeEdgeMutations();
  compare_treedags(treedag.View(), sample_tree.View());

  Merge<MADAG> merge{treedag.View().GetReferenceSequence()};
  merge.AddDAG(treedag.View());
  merge.AddDAG(sample_tree.View());
  merge.ComputeResultEdgeMutations();

  compare_treedags(treedag.View(), merge.GetResult());
}

[[maybe_unused]] static const auto test_added_write =
    add_test({test_write_protobuf, "Write protobuf"});

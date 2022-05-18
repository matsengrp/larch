#include "merge.hpp"

#include <iostream>
#include <fstream>

#include "test_common.hpp"
#include "history_dag_loader.hpp"
#include "benchmark.hpp"

void GetLabels(const HistoryDAG& tree, std::vector<NodeLabel>& labels,
               const std::string& refseq, const std::vector<CompactGenome>& mutations) {
  labels.resize(tree.GetNodes().size());
  for (auto iter : tree.TraversePreOrder()) {
    if (iter.IsRoot()) {
      continue;
    }
    const CompactGenome& muts = mutations.at(iter.GetEdge().GetId().value);
    NodeLabel& label = labels.at(iter.GetNode().GetId().value);
    const CompactGenome& parent_cgs =
        labels.at(iter.GetEdge().GetParent().GetId().value).first;
    label.first = parent_cgs;
    for (auto [pos, base] : muts) {
      if (base != refseq.at(pos - 1)) {
        label.first[pos] = base;
      } else {
        label.first.erase(pos);
      }
    }
  }

  for (Node node : tree.TraversePostOrder()) {
    if (node.IsLeaf()) {
      continue;
    }
    LeafSet& leaf_set = labels.at(node.GetId().value).second;
    for (auto clade : node.GetClades()) {
      std::set<CompactGenome> clade_leafs;
      for (Node child : clade | ranges::views::transform(Transform::GetChild)) {
        if (child.IsLeaf()) {
          clade_leafs.insert(labels.at(child.GetId().value).first);
        } else {
          for (auto& child_leafs : labels.at(child.GetId().value).second) {
            clade_leafs.insert(child_leafs.begin(), child_leafs.end());
          }
        }
      }
      leaf_set.insert(clade_leafs);
    }
  }
}

[[maybe_unused]] std::string ToString(const NodeLabel& label) {
  std::string result;
  for (auto [pos, mut] : label.first) {
    result += std::to_string(pos);
    result += mut;
    result += " ";
  }
  result += "[";
  for (auto& clade : label.second) {
    result += "(";
    for (auto& leaf : clade) {
      for (auto [pos, mut] : leaf) {
        result += std::to_string(pos);
        result += mut;
        result += " ";
      }
    }
    result += ") ";
  }
  result += "]";
  return result;
}

std::string ToDOT(Node node, const NodeLabel& label) {
  std::string result;
  size_t count = 0;
  for (auto [pos, mut] : label.first) {
    result += std::to_string(pos);
    result += mut;
    result += ++count % 3 == 0 ? "\\n" : " ";
  }
  result += "\\n[";
  result += std::to_string(node.GetId().value);
  result += "]";
  return result;
}

std::string ToDOT(const HistoryDAG& tree, const std::vector<NodeLabel>& labels) {
  std::string result;
  result += "digraph {\n";
  for (auto i : tree.GetEdges()) {
    std::string parent = ToDOT(i.GetParent(), labels.at(i.GetParent().GetId().value));
    std::string child = ToDOT(i.GetChild(), labels.at(i.GetChild().GetId().value));
    result += "  \"";
    result += parent;
    result += "\" -> \"";
    result += child;
    result += "\"\n";
  }
  result += "}";
  return result;
}

static void test_protobuf(const std::string& correct_path,
                          const std::vector<std::string>& paths) {
  std::vector<std::vector<CompactGenome>> mutations;
  std::vector<HistoryDAG> trees;
  for (auto& path : paths) {
    std::vector<CompactGenome> tree_mutations;
    trees.emplace_back(LoadHistoryDAGFromProtobufGZ(path, tree_mutations));
    mutations.emplace_back(std::move(tree_mutations));
  }
  std::string refseq;
  HistoryDAG correct_result = LoadHistoryDAGFromJsonGZ(correct_path, refseq);

  std::vector<std::vector<NodeLabel>> labels;
  std::vector<std::reference_wrapper<const HistoryDAG>> tree_refs;
  for (size_t i = 0; i < trees.size(); ++i) {
    std::vector<NodeLabel> tree_labels;
    GetLabels(trees.at(i), tree_labels, refseq, mutations.at(i));
    labels.emplace_back(std::move(tree_labels));
    tree_refs.push_back(trees.at(i));
  }

  Merge merged(refseq, std::move(tree_refs), labels);
  merged.Run();

  assert_equal(correct_result.GetNodes().size(), merged.GetResult().GetNodes().size(),
               "Nodes count");

  assert_equal(correct_result.GetEdges().size(), merged.GetResult().GetEdges().size(),
               "Edges count");
}

static void test_five_trees() {
  test_protobuf("data/test_5_trees/full_dag.json.gz",
                {
                    "data/test_5_trees/tree_0.pb.gz",
                    "data/test_5_trees/tree_1.pb.gz",
                    "data/test_5_trees/tree_2.pb.gz",
                    "data/test_5_trees/tree_3.pb.gz",
                    "data/test_5_trees/tree_4.pb.gz",
                });
}

static void test_case_2() {
  test_protobuf("data/testcase2/full_dag.json.gz", {
                                                       "data/testcase2/tree_0.pb.gz",
                                                       "data/testcase2/tree_1.pb.gz",
                                                       "data/testcase2/tree_2.pb.gz",
                                                       "data/testcase2/tree_3.pb.gz",
                                                       "data/testcase2/tree_4.pb.gz",
                                                   });
}

static void test_case_ref() {
  test_protobuf("data/testcaseref/full_dag.json.gz",
                {
                    "data/testcaseref/tree_0.pb.gz",
                    "data/testcaseref/tree_1.pb.gz",
                    "data/testcaseref/tree_2.pb.gz",
                    "data/testcaseref/tree_3.pb.gz",
                    "data/testcaseref/tree_4.pb.gz",
                });
}

[[maybe_unused]] static const auto test0_added = add_test({test_case_2, "Test case 2"});

[[maybe_unused]] static const auto test1_added =
    add_test({test_five_trees, "Merge 5 trees"});

[[maybe_unused]] static const auto test2_added =
    add_test({test_case_ref, "Tree with different ref"});

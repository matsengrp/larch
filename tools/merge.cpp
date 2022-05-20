#include <cstdlib>
#include <iostream>
#include <vector>

#include <range/v3/action/push_back.hpp>

#include "arguments.hpp"
#include "merge.hpp"
#include "history_dag_loader.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "merge [-z,--gzip] -i,--input file1 file2 ... [-o,--output filename]\n";
  std::cout << "  -i,--input     List of input files\n";
  std::cout << "  -z,--gzip      Input files are gzip compressed\n";
  std::cout << "  -o,--output    Save the output to filename (default is merged.pb)\n";

  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

static void GetLabels(const HistoryDAG& tree, std::vector<NodeLabel>& labels,
                      const std::string_view& refseq,
                      const std::vector<CompactGenome>& mutations) {
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

static HistoryDAG MergeTrees(const std::vector<std::string_view>& paths) {
  std::vector<std::vector<CompactGenome>> mutations;
  std::vector<HistoryDAG> trees;
  std::vector<std::string> ref_sequences;
  for (auto path : paths) {
    std::vector<CompactGenome> tree_mutations;
    std::string ref_seq;
    trees.emplace_back(LoadHistoryDAGFromProtobufGZ(path, ref_seq, tree_mutations));
    mutations.emplace_back(std::move(tree_mutations));
    ref_sequences.emplace_back(std::move(ref_seq));
  }

  std::vector<std::vector<NodeLabel>> labels;
  std::vector<std::reference_wrapper<const HistoryDAG>> tree_refs;
  for (size_t i = 0; i < trees.size(); ++i) {
    std::vector<NodeLabel> tree_labels;
    GetLabels(trees.at(i), tree_labels, ref_sequences.at(i), mutations.at(i));
    labels.emplace_back(std::move(tree_labels));
    tree_refs.push_back(trees.at(i));
  }

  Merge merge(ref_sequences.front(), std::move(tree_refs), labels);
  merge.Run();
  return merge.GetResult();
}

int main(int argc, char** argv) {
  Arguments args = GetArguments(argc, argv);

  bool gzip = false;
  std::vector<std::string_view> input_filenames;
  std::string result_filename = "merged.pb";

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "-i" or name == "--input") {
      ranges::action::push_back(input_filenames, params);
    } else if (name == "-z" or name == "--gzip") {
      gzip = true;
    } else if (name == "-o" or name == "--output") {
      if (params.empty()) {
        std::cerr << "Specify result file name.\n";
        Fail();
      }
      result_filename = *params.begin();
    }
  }

  if (input_filenames.size() < 2) {
    std::cerr << "Specify at least two input file names.\n";
    Fail();
  }

  std::cout << "Merging ";
  if (gzip) {
    std::cout << "gzip compressed ";
  }
  for (auto& i : input_filenames) {
    std::cout << i << " ";
  }
  std::cout << "into " << result_filename << "\n";

  HistoryDAG result = MergeTrees(input_filenames);

  std::cout << "DAG nodes: " << result.GetNodes().size() << "\n";
  std::cout << "DAG edges: " << result.GetEdges().size() << "\n";

  return EXIT_SUCCESS;
}

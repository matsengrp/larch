#include "larch/merge/merge.hpp"

#include <iostream>
#include <filesystem>
#include <vector>

#include "test_common.hpp"
#include "larch/dag_loader.hpp"

static void test_linearham_load_and_merge() {
  const std::string data_dir = "data/linearham/";
  const std::vector<std::string_view> fasta_paths = {
      "data/linearham/igh.fa", "data/linearham/naive-seqs.fa"};
  const std::string ref_path = data_dir + "reference_sequence.txt";

  // Collect all -rerooted.treefile paths
  std::vector<std::string> tree_paths;
  for (const auto& entry : std::filesystem::directory_iterator{data_dir}) {
    if (entry.path().extension() != ".treefile") {
      continue;
    }
    auto filename = entry.path().filename().string();
    if (filename.find("-rerooted") == std::string::npos) {
      continue;
    }
    tree_paths.push_back(entry.path().string());
  }
  std::sort(tree_paths.begin(), tree_paths.end());

  TestAssert(tree_paths.size() > 0);
  std::cout << "Found " << tree_paths.size() << " tree files." << std::endl;

  // Load reference sequence for Merge
  std::string ref_seq = LoadReferenceSequence(ref_path);
  TestAssert(!ref_seq.empty());

  // Load all trees
  std::vector<MADAGStorage<>> trees;
  trees.reserve(tree_paths.size());
  for (const auto& tree_path : tree_paths) {
    trees.push_back(LoadTreeFromFastaNewick(fasta_paths, tree_path, ref_path));
    auto view = trees.back().View();
    // Verify tree was loaded correctly
    TestAssert(view.GetNodesCount() > 0);
    TestAssert(view.GetEdgesCount() > 0);
    // Check that leaf nodes have SampleIds
    for (auto node : view.GetNodes()) {
      if (node.IsLeaf()) {
        TestAssert(node.HaveSampleId());
      }
    }
    // Compute compact genomes (required for merge)
    view.RecomputeCompactGenomes(true);
  }

  std::cout << "Loaded " << trees.size() << " trees." << std::endl;

  // Create views for merge
  std::vector<MADAG> tree_views;
  tree_views.reserve(trees.size());
  for (auto& tree : trees) {
    tree_views.push_back(tree.View());
  }

  // Merge all trees
  Merge merge(ref_seq);
  merge.AddDAGs(tree_views);
  merge.GetResult().GetRoot().Validate(true, true);

  MADAGToDOT(merge.GetResult(), std::cout);

  // Basic sanity checks on merged result
  TestAssert(merge.GetResult().GetNodesCount() > 0);
  TestAssert(merge.GetResult().GetEdgesCount() > 0);

  std::cout << "Merged DAG has " << merge.GetResult().GetNodesCount() << " nodes and "
            << merge.GetResult().GetEdgesCount() << " edges." << std::endl;
}

[[maybe_unused]] static const auto test_added =
    add_test({test_linearham_load_and_merge, "LinearHam: Load and merge trees"});

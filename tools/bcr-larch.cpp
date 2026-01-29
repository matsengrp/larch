#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <vector>

#ifdef DISABLE_PARALLELISM
#include <tbb/global_control.h>
#endif

#include "larch/merge/merge.hpp"
#include "larch/dag_loader.hpp"
#include "tools_common.hpp"

/*
Usage:
  bcr-larch --fasta FILE --fasta FILE ... --trees DIR --reference FILE -o FILE

  Options:
  - --fasta FILE - Path to FASTA file (can be specified multiple times)
  - --trees DIR - Directory containing Newick tree files
  - --reference FILE - Path to reference sequence file
  - -o,--output FILE - Path to output DAG protobuf file
  - --tree-suffix SUFFIX - Suffix for tree files (default: -rerooted.treefile)

  Example:
  ./bin/bcr-larch \
    --fasta data/linearham/igh.fa \
    --fasta data/linearham/naive-seqs.fa \
    --trees data/linearham \
    --reference data/linearham/reference_sequence.txt \
    -o output.pb

  The tool:
  1. Loads sequences from multiple FASTA files
  2. Finds tree files matching the suffix in the specified directory
  3. Loads each tree using LoadTreeFromFastaNewick
  4. Merges all trees using the Merge class
  5. Saves the result DAG in protobuf format
*/

[[noreturn]] static void Usage() {
  std::string program_desc =
      "bcr-larch: Load and merge trees from FASTA sequences and Newick topologies";

  std::vector<std::string> usage_examples = {
      {"bcr-larch --fasta FILE --fasta FILE ... --trees DIR --reference FILE -o FILE"}};

  std::vector<std::pair<std::string, std::string>> flag_desc_pairs = {
      {"--fasta FILE",
       "Path to FASTA file (can be specified multiple times) (REQUIRED)\n"
       "Typically: leaf sequences and naive/root sequences"},
      {"--trees DIR", "Directory containing Newick tree files (*.treefile) (REQUIRED)"},
      {"--reference FILE", "Path to reference sequence file (REQUIRED)"},
      {"-o,--output FILE", "Path to output DAG protobuf file (REQUIRED)"},
      {"--tree-suffix SUFFIX",
       "Suffix for tree files to load (default: -rerooted.treefile)"},
  };

  std::cout << FormatUsage(program_desc, usage_examples, flag_desc_pairs);

  std::exit(EXIT_SUCCESS);
}

int main(int argc, char** argv) try {
#ifdef DISABLE_PARALLELISM
  tbb::global_control gc(tbb::global_control::max_allowed_parallelism, 1);
#endif
  Arguments args = GetArguments(argc, argv);

  std::vector<std::string> fasta_paths;
  std::string trees_dir;
  std::string reference_path;
  std::string output_path;
  std::string tree_suffix = "-rerooted.treefile";

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "--version") {
      Version("bcr-larch");
    } else if (name == "--fasta") {
      std::string path;
      ParseOption(name, params, path, 1);
      fasta_paths.push_back(path);
    } else if (name == "--trees") {
      ParseOption(name, params, trees_dir, 1);
    } else if (name == "--reference") {
      ParseOption(name, params, reference_path, 1);
    } else if (name == "-o" or name == "--output") {
      ParseOption(name, params, output_path, 1);
    } else if (name == "--tree-suffix") {
      ParseOption(name, params, tree_suffix, 1);
    } else {
      std::cerr << "Unknown argument '" << name << "'.\n";
      Fail();
    }
  }

  // Validate required arguments
  if (fasta_paths.empty()) {
    std::cerr << "ERROR: At least one FASTA file is required (--fasta).\n";
    Fail();
  }
  if (trees_dir.empty()) {
    std::cerr << "ERROR: Trees directory is required (--trees).\n";
    Fail();
  }
  if (reference_path.empty()) {
    std::cerr << "ERROR: Reference sequence file is required (--reference).\n";
    Fail();
  }
  if (output_path.empty()) {
    std::cerr << "ERROR: Output file is required (-o/--output).\n";
    Fail();
  }

  // Convert to string_view for LoadTreeFromFastaNewick
  std::vector<std::string_view> fasta_paths_sv;
  fasta_paths_sv.reserve(fasta_paths.size());
  for (const auto& path : fasta_paths) {
    fasta_paths_sv.push_back(path);
  }

  // Collect tree files from directory
  std::vector<std::string> tree_paths;
  for (const auto& entry : std::filesystem::directory_iterator{trees_dir}) {
    if (entry.path().extension() != ".treefile") {
      continue;
    }
    auto filename = entry.path().filename().string();
    if (filename.find(tree_suffix.substr(0, tree_suffix.find('.'))) ==
        std::string::npos) {
      continue;
    }
    tree_paths.push_back(entry.path().string());
  }
  std::sort(tree_paths.begin(), tree_paths.end());

  if (tree_paths.empty()) {
    std::cerr << "ERROR: No tree files found in '" << trees_dir << "' with suffix '"
              << tree_suffix << "'.\n";
    Fail();
  }

  std::cout << "Found " << tree_paths.size() << " tree files.\n";
  std::cout << "Loading FASTA files:";
  for (const auto& path : fasta_paths) {
    std::cout << " " << path;
  }
  std::cout << "\n";

  // Load reference sequence
  std::string ref_seq = LoadReferenceSequence(reference_path);
  if (ref_seq.empty()) {
    std::cerr << "ERROR: Failed to load reference sequence from '" << reference_path
              << "'.\n";
    Fail();
  }
  std::cout << "Reference sequence: " << ref_seq.size() << " bases.\n";

  // Load all trees
  std::vector<MADAGStorage<>> trees;
  trees.reserve(tree_paths.size());
  for (const auto& tree_path : tree_paths) {
    std::cout << "Loading: " << tree_path << "\n";
    trees.push_back(LoadTreeFromFastaNewick(fasta_paths_sv, tree_path, reference_path));
    auto view = trees.back().View();
    view.RecomputeCompactGenomes(true);
  }
  std::cout << "Loaded " << trees.size() << " trees.\n";

  // Create views for merge
  std::vector<MADAG> tree_views;
  tree_views.reserve(trees.size());
  for (auto& tree : trees) {
    tree_views.push_back(tree.View());
  }

  // Merge all trees
  std::cout << "Merging trees...\n";
  Merge merge(ref_seq);
  merge.AddDAGs(tree_views);
  merge.ComputeResultEdgeMutations();
  merge.GetResult().GetRoot().Validate(true, true);

  std::cout << "Merged DAG has " << merge.GetResult().GetNodesCount() << " nodes and "
            << merge.GetResult().GetEdgesCount() << " edges.\n";

  // Save to protobuf
  std::cout << "Saving to: " << output_path << "\n";
  StoreDAG(merge.GetResult(), output_path, FileFormat::ProtobufDAG);

  std::cout << "Done.\n";
  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << "Uncaught exception: " << e.what() << std::endl;
  std::terminate();
} catch (...) {
  std::abort();
}

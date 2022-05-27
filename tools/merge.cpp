#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>

#include <range/v3/action/push_back.hpp>

#include "arguments.hpp"
#include "merge.hpp"
#include "history_dag_loader.hpp"
#include "benchmark.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "merge [-r,--refseq file] -i,--input file1 file2 ... "
               "[-o,--output filename]\n";
  std::cout << "  -i,--input     List of input files\n";
  std::cout << "  -o,--output    Save the output to filename (default is merged.pb)\n";
  std::cout << "  -r,--refseq    Read reference sequence from Json file\n";

  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

static int MergeTrees(const std::vector<std::string_view>& paths,
                      std::string_view refseq_json_path, std::string_view out_path) {
  std::vector<std::vector<Mutations>> mutations;
  std::vector<HistoryDAG> trees;
  std::string reference_sequence;

  reference_sequence = LoadRefseqFromJsonGZ(refseq_json_path);

  trees.resize(paths.size());
  mutations.resize(paths.size());
  std::vector<std::pair<size_t, std::string_view>> paths_idx;
  for (size_t i = 0; i < paths.size(); ++i) {
    paths_idx.push_back({i, paths.at(i)});
  }
  std::cout << "Loading trees ";
  tbb::parallel_for_each(paths_idx.begin(), paths_idx.end(), [&](auto path_idx) {
    std::vector<Mutations> tree_mutations;
    std::cout << "." << std::flush;
    trees.at(path_idx.first) = LoadTreeFromProtobufGZ(path_idx.second, tree_mutations);
    mutations.at(path_idx.first) = std::move(tree_mutations);
  });
  std::cout << " done."
            << "\n";

  Benchmark merge_time;
  Merge merge(reference_sequence, trees, mutations);
  merge_time.start();
  merge.Run();
  merge_time.stop();
  std::cout << "\nDAGs merged in " << merge_time.durationMs() << " ms\n";

  std::cout << "DAG nodes: " << merge.GetResult().GetNodes().size() << "\n";
  std::cout << "DAG edges: " << merge.GetResult().GetEdges().size() << "\n";

  StoreDAGToProtobuf(merge.GetResult(), reference_sequence,
                     merge.CalculateResultEdgeMutations(), out_path);

  return EXIT_SUCCESS;
}

int main(int argc, char** argv) {
  Arguments args = GetArguments(argc, argv);

  std::vector<std::string_view> input_filenames;
  std::string result_filename = "merged.pb";
  std::string refseq_filename;

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "-i" or name == "--input") {
      ranges::action::push_back(input_filenames, params);
    } else if (name == "-o" or name == "--output") {
      if (params.empty()) {
        std::cerr << "Specify result file name.\n";
        Fail();
      }
      result_filename = *params.begin();
    } else if (name == "-r" or name == "--refseq") {
      if (params.empty()) {
        std::cerr << "Specify reference sequence Json file name.\n";
        Fail();
      }
      refseq_filename = *params.begin();
    }
  }

  if (input_filenames.size() < 2) {
    std::cerr << "Specify at least two input file names.\n";
    Fail();
  }

  return MergeTrees(input_filenames, refseq_filename, result_filename);
}

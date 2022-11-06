#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>

#include "arguments.hpp"
#include "larch/merge/merge.hpp"
#include "larch/dag_loader.hpp"
#include "benchmark.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "merge [-r,--refseq file] [-d,--dag] -i,--input file1 file2 ... "
               "[-o,--output filename]\n";
  std::cout << "  -i,--input     List of input files\n";
  std::cout << "  -o,--output    Save the output to filename (default is merged.pb)\n";
  std::cout << "  -r,--refseq    Read reference sequence from Json file\n";
  std::cout << "  -d,--dag       Input files are DAGs\n";

  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

static int MergeTrees(const std::vector<std::string_view>& paths,
                      std::string_view refseq_json_path, std::string_view out_path,
                      bool dags) {
  std::vector<MADAGStorage> trees;
  std::string reference_sequence =
      std::string{LoadDAGFromJson(refseq_json_path).View().GetReferenceSequence()};

  trees.resize(paths.size());
  std::vector<std::pair<size_t, std::string_view>> paths_idx;
  for (size_t i = 0; i < paths.size(); ++i) {
    paths_idx.emplace_back(i, paths.at(i));
  }
  std::cout << "Loading trees ";
  tbb::parallel_for_each(paths_idx.begin(), paths_idx.end(), [&](auto path_idx) {
    std::cout << "." << std::flush;
    trees.at(path_idx.first) =
        dags ? LoadDAGFromProtobuf(path_idx.second)
             : LoadTreeFromProtobuf(path_idx.second, reference_sequence);
  });
  std::cout << " done."
            << "\n";

  Benchmark merge_time;
  Merge merge(reference_sequence);
  std::vector<MADAG> tree_refs{trees.begin(), trees.end()};
  merge_time.start();
  merge.AddDAGs(tree_refs);
  merge_time.stop();
  std::cout << "\nDAGs merged in " << merge_time.durationMs() << " ms\n";

  std::cout << "DAG nodes: " << merge.GetResult().GetNodesCount() << "\n";
  std::cout << "DAG edges: " << merge.GetResult().GetEdgesCount() << "\n";

  StoreDAGToProtobuf(merge.GetResult(), out_path);

  return EXIT_SUCCESS;
}

int main(int argc, char** argv) {
  Arguments args = GetArguments(argc, argv);

  std::vector<std::string_view> input_filenames;
  std::string result_filename = "merged.pb";
  std::string refseq_filename;
  bool dags = false;

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "-i" or name == "--input") {
      ranges::actions::push_back(input_filenames, params);
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
    } else if (name == "-d" or name == "--dag") {
      dags = true;
    }
  }

  if (input_filenames.size() < 2) {
    std::cerr << "Specify at least two input file names.\n";
    Fail();
  }

  return MergeTrees(input_filenames, refseq_filename, result_filename, dags);
}

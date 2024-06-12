#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>

#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score_binary.hpp"
#include "arguments.hpp"
#include "larch/rf_distance.hpp"
#include "larch/merge/merge.hpp"
#include "larch/dag_loader.hpp"
#include "larch/benchmark.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "merge [-r,--refseq file] -i,--input infile1 infile2 ... "
               "[-o,--output outfile]\n";
  std::cout << "  -i,--input       Paths to input Trees/DAGs\n";
  std::cout << "  -o,--output      Path to output DAG (default: merged.pb)\n";
  std::cout << "  -r,--refseq      Read reference sequence from Json file\n";
  std::cout << "  -t,--trim        Trim output (default: best parsimony)\n";
  std::cout << "  --rf             Trim output to minimize RF distance to provided "
               "DAG file\n";
  std::cout << "  -s,--sample      Sample a single tree from DAG\n";
  std::cout << "  --input-format   List the input file formats (default: inferred)\n";
  std::cout << "  --output-format  Output file format (default: inferred)\n";
  std::cout << "  --rf-format      RF file format (default: inferred)\n";
  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

static void MergeTrees(const std::vector<std::string_view>& input_paths,
                       const std::vector<FileFormat>& input_formats,
                       std::string refseq_path, std::string_view output_path,
                       FileFormat output_format, bool trim, bool sample_tree,
                       std::string rf_path, FileFormat rf_format) {
  std::vector<MADAGStorage<>> trees;
  std::vector<size_t> trees_id;

  trees.reserve(input_paths.size());
  for (size_t i = 0; i < input_paths.size(); ++i) {
    trees_id.push_back(i);
    trees.push_back(MADAGStorage<>::EmptyDefault());
  }
  std::cout << "Loading trees ";
  ParallelForEach(trees_id, [&](auto tree_id) {
    const auto tree_path = input_paths.at(tree_id);
    const auto tree_format = input_formats.at(tree_id);
    std::cout << "." << std::flush;
    trees.at(tree_id) = LoadDAG(tree_path, tree_format, refseq_path);
    trees.at(tree_id).View().RecomputeCompactGenomes();
  });
  std::cout << " done.\n";

  Benchmark merge_time;
  Merge merge(trees.front().View().GetReferenceSequence());
  std::vector<MADAG> tree_refs{trees.begin(), trees.end()};
  merge_time.start();
  merge.AddDAGs(tree_refs);
  merge_time.stop();
  std::cout << "\nDAGs merged in " << merge_time.durationMs() << " ms\n";

  std::cout << "DAG leave(without trimming): " << merge.GetResult().GetLeafsCount()
            << "\n";
  std::cout << "DAG nodes(without trimming): " << merge.GetResult().GetNodesCount()
            << "\n";
  std::cout << "DAG edges(without trimming): " << merge.GetResult().GetEdgesCount()
            << "\n";

  merge.ComputeResultEdgeMutations();
  if (trim) {
    if (rf_path.empty()) {
      SubtreeWeight<BinaryParsimonyScore, MergeDAG> weight{merge.GetResult()};
      if (sample_tree) {
        std::cout << "sampling a tree from the minweight options\n";
        StoreDAG(weight.MinWeightSampleTree({}).View(), output_path, output_format);
      } else {
        StoreDAG(weight.TrimToMinWeight({}).View(), output_path, output_format);
      }
    } else {
      auto tree = LoadDAG(rf_path, rf_format, refseq_path);
      Merge comparetree(tree.View().GetReferenceSequence());
      comparetree.AddDAG(tree.View());
      comparetree.ComputeResultEdgeMutations();
      SubtreeWeight<SumRFDistance, MergeDAG> min_sum_rf_dist{merge.GetResult()};
      SumRFDistance min_rf_weight_ops{comparetree, merge};
      if (sample_tree) {
        std::cout << "sampling a tree from the minweight options\n";
        StoreDAG(min_sum_rf_dist.MinWeightSampleTree(min_rf_weight_ops, {}).View(),
                 output_path, output_format);
      } else {
        StoreDAG(min_sum_rf_dist.TrimToMinWeight(min_rf_weight_ops).View(), output_path,
                 output_format);
      }
    }
  } else {
    if (sample_tree) {
      std::cout << "sampling a tree from the merge DAG\n";
      SubtreeWeight<BinaryParsimonyScore, MergeDAG> weight{merge.GetResult()};
      StoreDAG(weight.SampleTree({}).View(), output_path, output_format);
    } else {
      StoreDAG(merge.GetResult(), output_path, output_format);
    }
  }
}

int main(int argc, char** argv) try {
  Arguments args = GetArguments(argc, argv);

  std::vector<std::string_view> input_paths;
  std::vector<FileFormat> input_formats;
  std::string output_path = "merged.dagbin";
  FileFormat output_format = FileFormat::Infer;
  std::string refseq_path;
  std::string rf_path;
  FileFormat rf_format = FileFormat::Infer;
  bool trim = false;
  bool sample_tree = false;

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "-i" or name == "--input") {
      ranges::actions::push_back(input_paths, params);
    } else if (name == "-o" or name == "--output") {
      if (params.empty()) {
        std::cerr << "Specify result file name.\n";
        Fail();
      }
      output_path = *params.begin();
    } else if (name == "-r" or name == "--refseq") {
      if (params.empty()) {
        std::cerr << "Specify reference sequence Json file name.\n";
        Fail();
      }
      refseq_path = *params.begin();
    } else if (name == "-t" or name == "--trim") {
      trim = true;
    } else if (name == "--rf") {
      if (params.empty()) {
        std::cerr << "Specify rf-trim protobuf file name.\n";
        Fail();
      }
      rf_path = *params.begin();
    } else if (name == "-s" or name == "--sample") {
      sample_tree = true;
    } else if (name == "--input-format") {
      if (params.empty()) {
        std::cerr << "Specify input formats.\n";
        Fail();
      }
      for (auto param : params) {
        input_formats.push_back(InferFileFormat(param));
      }
    } else if (name == "--output-format") {
      if (params.empty()) {
        std::cerr << "Specify output format.\n";
        Fail();
      }
      output_format = InferFileFormat(*params.begin());
    } else if (name == "--rf-format") {
      if (params.empty()) {
        std::cerr << "Specify output format.\n";
        Fail();
      }
      rf_format = InferFileFormat(*params.begin());
    } else {
      std::cerr << "Unknown argument '" << name << "'.\n";
      Fail();
    }
  }

  if (input_paths.size() < 1) {
    std::cerr << "Specify at least one input file names.\n";
    Fail();
  }
  for (auto input_path : input_paths) {
    std::cout << input_path << " to be merged\n";
  }
  if (input_formats.empty()) {
    for (int i = 0; i < int(input_paths.size()); i++) {
      input_formats.push_back(FileFormat::Infer);
    }
  }
  if (input_paths.size() != input_formats.size()) {
    std::cerr << "Specify format for each input.\n";
    Fail();
  }

  bool is_input_dag = refseq_path.empty();
  for (size_t i = 0; i < input_formats.size(); i++) {
    auto& input_format = input_formats.at(i);
    const auto& input_path = input_paths.at(i);
    if (input_format == FileFormat::Infer) {
      input_format = InferFileFormat(input_path);
      if (input_format == FileFormat::Protobuf) {
        input_format =
            is_input_dag ? FileFormat::ProtobufTree : FileFormat::ProtobufDAG;
      }
    }
  }
  if (output_format == FileFormat::Infer) {
    output_format = InferFileFormat(output_path);
    if (output_format == FileFormat::Protobuf) {
      output_format = FileFormat::ProtobufDAG;
    }
  }
  if (!rf_path.empty() and rf_format == FileFormat::Infer) {
    rf_format = InferFileFormat(rf_path);
    if (rf_format == FileFormat::Protobuf) {
      rf_format = FileFormat::ProtobufDAG;
    }
  }

  MergeTrees(input_paths, input_formats, refseq_path, output_path, output_format, trim,
             sample_tree, rf_path, rf_format);
  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << "Uncaught exception: " << e.what() << std::endl;
  std::terminate();
} catch (...) {
  std::abort();
}

#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>

#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/weight_accumulator.hpp"
#include "larch/subtree/parsimony_score_binary.hpp"
#include "tools_common.hpp"
#include "larch/rf_distance.hpp"
#include "larch/merge/merge.hpp"
#include "larch/dag_loader.hpp"
#include "larch/benchmark.hpp"

[[noreturn]] static void Usage() {
  std::string program_desc =
      "dag-util: General utility for manipulating (e.g. combining, pruning) or "
      "inspecting DAGs/trees";

  std::vector<std::string> usage_examples = {
      {"dag-util [-r,--refseq FILE] -i,--input FILE1 FILE2 ... [-o,--output FILE]"}};

  std::vector<std::pair<std::string, std::string>> flag_desc_pairs = {
      {"-i,--input FILE [...]", "Paths to input DAG/Tree files (REQUIRED)"},
      {"-o,--output FILE",
       "Path to output DAG file (default: does not save result DAG)"},
      {"-r,--MAT-refseq-file FILE",
       "Path to json reference sequence file \n"
       "(REQUIRED if input file is a MAT protobuf)"},
      {"-t,--trim", "Trim output (default: best parsimony)"},
      {"--rf FILE", "Trim output to minimize RF distance to provided DAG file"},
      {"-s,--sample", "Sample a single tree from DAG"},
      {"--dag-info", "Print DAG info (parsimony scores, sum RF distances)"},
      {"--parsimony", "Print all DAG parsimony scores"},
      {"--sum-rf-distance", "Print all DAG sum RF distances"},
      {"--input-format ENUM [...]",
       "Specify input file formats (default: inferred) \n"
       "[dagbin, dag-pb, tree-pb, dag-json]"},
      {"--output-format ENUM",
       "Specify output file format (default: inferred) \n"
       "[dagbin, dag-pb]"},
      {"--rf-format ENUM",
       "Specify RF file format (default: inferred) \n"
       "[dagbin, dag-pb, tree-pb, dag-json]"},
  };

  std::cout << FormatUsage(program_desc, usage_examples, flag_desc_pairs);

  std::exit(EXIT_SUCCESS);
}

static void MergeTrees(const std::vector<std::string_view>& input_paths,
                       const std::vector<FileFormat>& input_formats,
                       std::string refseq_path, std::string_view output_path,
                       FileFormat output_format, bool trim, bool sample_tree,
                       std::string rf_path, FileFormat rf_format,
                       bool do_print_dag_info, bool do_print_parsimony,
                       bool do_print_rf_distance) {
  std::vector<MADAGStorage<>> trees;
  std::vector<size_t> trees_id;

  trees.reserve(input_paths.size());
  for (size_t i = 0; i < input_paths.size(); ++i) {
    trees_id.push_back(i);
    trees.push_back(MADAGStorage<>::EmptyDefault());
  }
  std::cout << "Loading trees... ";
  ParallelForEach(trees_id, [&](auto tree_id) {
    const auto tree_path = input_paths.at(tree_id);
    const auto tree_format = input_formats.at(tree_id);
    std::cout << " . " << std::flush;
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

  if (do_print_dag_info) {
    auto scorecount_compare = [](const auto& lhs, const auto& rhs) {
      return lhs.first < rhs.first;
    };
    auto root_node = merge.GetResult().GetRoot();

    SubtreeWeight<TreeCount, MergeDAG> tree_counter{merge.GetResult()};
    auto tree_count = tree_counter.ComputeWeightBelow(root_node, {});
    std::cout << "tree_count: " << tree_count << std::endl;

    if (do_print_parsimony) {
      // Compute all Parsimony Scores
      SubtreeWeight<WeightAccumulator<BinaryParsimonyScore>, MergeDAG> parsimony_scorer{
          merge.GetResult()};
      auto all_parsimony_results = parsimony_scorer.ComputeWeightBelow(root_node, {});
      const auto& all_parsimony_scores = all_parsimony_results.GetWeights();

      auto min_parsimony_score = *std::min_element(
          all_parsimony_scores.begin(), all_parsimony_scores.end(), scorecount_compare);
      auto max_parsimony_score = *std::max_element(
          all_parsimony_scores.begin(), all_parsimony_scores.end(), scorecount_compare);

      std::cout << "parsimony_all: " << all_parsimony_scores.size() << "\n"
                << all_parsimony_scores << "\n";
      std::cout << "parsimony_min: score:" << min_parsimony_score.first
                << ", count:" << min_parsimony_score.second << "\n";
      std::cout << "parsimony_max: score:" << max_parsimony_score.first
                << ", count:" << max_parsimony_score.second << " \n";
    }

    if (do_print_rf_distance) {
      // Compute all Sum RF Distances
      SubtreeWeight<WeightAccumulator<SumRFDistance>, MergeDAG> sum_rf_dist_scorer{
          merge.GetResult()};
      SumRFDistance sum_rf_dist_weight_ops{merge, merge};
      auto all_sum_rf_dist_results =
          sum_rf_dist_scorer.ComputeWeightBelow(root_node, {sum_rf_dist_weight_ops});
      auto all_sum_rf_dist_scores = all_sum_rf_dist_results.GetWeights().Copy();
      auto shift_sum = sum_rf_dist_weight_ops.GetOps().GetShiftSum();
      for (auto& score_count : all_sum_rf_dist_scores) {
        score_count.first += shift_sum;
      }

      auto min_sum_rf_dist_scores =
          *std::min_element(all_sum_rf_dist_scores.begin(),
                            all_sum_rf_dist_scores.end(), scorecount_compare);
      auto max_sum_rf_dist_scores =
          *std::max_element(all_sum_rf_dist_scores.begin(),
                            all_sum_rf_dist_scores.end(), scorecount_compare);

      std::cout << "sum_rf_dist_all: " << all_sum_rf_dist_scores.size() << "\n"
                << all_sum_rf_dist_scores << "\n";
      std::cout << "sum_rf_dist_min: score:" << min_sum_rf_dist_scores.first
                << ", count:" << min_sum_rf_dist_scores.second << "\n";
      std::cout << "sum_rf_dist_max: score:" << max_sum_rf_dist_scores.first
                << ", count:" << max_sum_rf_dist_scores.second << "\n";
    }
  }

  if (!output_path.empty()) {
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
          StoreDAG(min_sum_rf_dist.TrimToMinWeight(min_rf_weight_ops).View(),
                   output_path, output_format);
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
}

int main(int argc, char** argv) try {
  Arguments args = GetArguments(argc, argv);

  std::vector<std::string_view> input_paths;
  std::vector<FileFormat> input_formats;
  std::string output_path;
  FileFormat output_format = FileFormat::Infer;
  std::string refseq_path;
  std::string rf_path;
  FileFormat rf_format = FileFormat::Infer;
  bool trim = false;
  bool sample_tree = false;
  bool do_print_dag_info = false;
  bool do_print_parsimony = false;
  bool do_print_rf_distance = false;

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "--version") {
      Version();
    } else if (name == "-i" or name == "--input") {
      ParseOption<false>(name, params, input_paths, -1);
      ranges::actions::push_back(input_paths, params);
    } else if (name == "-o" or name == "--output") {
      ParseOption(name, params, output_path, 1);
    } else if (name == "-r" or name == "--refseq") {
      ParseOption(name, params, refseq_path, 1);
    } else if (name == "-t" or name == "--trim") {
      ParseOption<false>(name, params, trim, 0);
      trim = true;
    } else if (name == "--rf") {
      ParseOption(name, params, rf_path, 1);
    } else if (name == "-s" or name == "--sample") {
      ParseOption<false>(name, params, sample_tree, 0);
      sample_tree = true;
    } else if (name == "--dag-info") {
      ParseOption<false>(name, params, do_print_dag_info, 0);
      do_print_dag_info = true;
      do_print_parsimony = true;
      do_print_rf_distance = true;
    } else if (name == "--print-pars") {
      ParseOption<false>(name, params, do_print_parsimony, 0);
      do_print_dag_info = true;
      do_print_parsimony = true;
    } else if (name == "--print-rf") {
      ParseOption<false>(name, params, do_print_rf_distance, 0);
      do_print_dag_info = true;
      do_print_rf_distance = true;
    } else if (name == "--input-format") {
      ParseOption<false>(name, params, input_formats, -1);
      for (auto param : params) {
        input_formats.push_back(InferFileFormat(param));
      }
    } else if (name == "--output-format") {
      std::string temp;
      ParseOption(name, params, temp, 1);
      output_format = InferFileFormat(temp);
    } else if (name == "--rf-format") {
      std::string temp;
      ParseOption(name, params, temp, 1);
      rf_format = InferFileFormat(temp);
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
    std::cerr << "Specify input format for each input file.\n";
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
            is_input_dag ? FileFormat::ProtobufDAG : FileFormat::ProtobufTree;
      }
    }
  }
  if (!output_path.empty() and output_format == FileFormat::Infer) {
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
             sample_tree, rf_path, rf_format, do_print_dag_info, do_print_parsimony,
             do_print_rf_distance);
  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << "Uncaught exception: " << e.what() << std::endl;
  std::terminate();
} catch (...) {
  std::abort();
}

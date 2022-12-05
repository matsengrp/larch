#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <chrono>

#include <unistd.h>
#include <sys/wait.h>

#include "arguments.hpp"
#include "larch/dag_loader.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/weight_accumulator.hpp"
#include "larch/subtree/tree_count.hpp"
#include "larch/subtree/parsimony_score_binary.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/merge/merge.hpp"
#include "benchmark.hpp"
#include <mpi.h>

#include "larch/usher_glue.hpp"

MADAGStorage optimize_dag_direct(MADAG dag, Move_Found_Callback& callback);
[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "larch-usher -i,--input file -o,--output file [-m,--matopt file] "
               "[-c,--count number]\n";
  std::cout << "  -i,--input   Path to input DAG\n";
  std::cout << "  -o,--output  Path to output DAG\n";
  std::cout << "  -m,--matopt  Path to matOptimize executable. Default: matOptimize\n";
  std::cout << "  -l,--logpath Path for logging\n";
  std::cout << "  -c,--count   Number of iterations. Default: 1\n";
  std::cout
      << "  --move-coeff-nodes   New node coefficient for scoring moves. Default: 1\n";
  std::cout << "  --move-coeff-pscore  Parsimony score coefficient for scoring moves. "
               "Default: 1\n";
  std::cout << "  --sample-best-tree   Only sample trees with best achieved parsimony "
               "score.\n";
  std::cout << "  -r,--MAT-refseq-file   Provide a path to a file containing a "
               "reference sequence\nif input points to MAT protobuf\n";

  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

static int ParseNumber(std::string_view str) {
  int result{};
  std::istringstream stream{std::string{str}};
  stream >> result;
  if (stream.fail()) {
    throw std::runtime_error("Invalid number");
  }
  return result;
}

void check_edge_mutations(MADAG madag);

std::vector<std::vector<const CompactGenome*>> clades_union(
    const std::vector<std::vector<const CompactGenome*>>& lhs,
    const std::vector<std::vector<const CompactGenome*>>& rhs) {
  std::vector<std::vector<const CompactGenome*>> result;

  for (auto [lhs_clade, rhs_clade] : ranges::views::zip(lhs, rhs)) {
    std::vector<const CompactGenome*> clade{lhs_clade};
    clade.insert(clade.end(), rhs_clade.begin(), rhs_clade.end());
    ranges::sort(clade);
    ranges::unique(clade);
    result.push_back(std::move(clade));
  }

  ranges::sort(result);
  return result;
}

std::vector<std::vector<const CompactGenome*>> clades_difference(
    const std::vector<std::vector<const CompactGenome*>>& lhs,
    const std::vector<std::vector<const CompactGenome*>>& rhs) {
  std::vector<std::vector<const CompactGenome*>> result;

  for (auto [lhs_clade, rhs_clade] : ranges::views::zip(lhs, rhs)) {
    std::vector<const CompactGenome*> clade;
    std::set_difference(lhs_clade.begin(), lhs_clade.end(), rhs_clade.begin(),
                        rhs_clade.end(), std::inserter(clade, clade.begin()));
    ranges::sort(clade);
    ranges::unique(clade);
    result.push_back(std::move(clade));
  }

  ranges::sort(result);
  return result;
}

// NOLINTNEXTLINE(cppcoreguidelines-virtual-class-destructor)
struct Larch_Move_Found_Callback : public Move_Found_Callback {
  Larch_Move_Found_Callback(const Merge& merge, MADAG sample,
                            const std::vector<NodeId>& sample_dag_ids)
      : merge_{merge},
        sample_{sample},
        sample_dag_ids_{sample_dag_ids},
        move_score_coeffs_{1, 1} {}
  Larch_Move_Found_Callback(
      const Merge& merge, MADAG sample, const std::vector<NodeId>& sample_dag_ids,
      std::pair<int, int> move_score_coeffs)  // NOLINT(modernize-pass-by-value)
      : merge_{merge},
        sample_{sample},
        sample_dag_ids_{sample_dag_ids},
        move_score_coeffs_{move_score_coeffs} {}
  bool operator()(Profitable_Moves& move, int /* best_score_change */,
                  std::vector<Node_With_Major_Allele_Set_Change>&
                  /* node_with_major_allele_set_change */) override {
    int new_nodes_count = 0;
    if (move_score_coeffs_.first != 0) {
      NodeId src_id = sample_dag_ids_.at(move.src->node_id);
      NodeId dst_id = sample_dag_ids_.at(move.dst->node_id);
      NodeId lca_id = sample_dag_ids_.at(move.LCA->node_id);

      const auto& src_clades =
          merge_.GetResultNodeLabels().at(src_id.value).GetLeafSet()->GetClades();
      const auto& dst_clades =
          merge_.GetResultNodeLabels().at(dst_id.value).GetLeafSet()->GetClades();

      MAT::Node* curr_node = move.src;
      while (not(curr_node->node_id == lca_id.value)) {
        MADAG::NodeView node = merge_.GetResult().Get(NodeId{curr_node->node_id});
        const auto& clades = merge_.GetResultNodeLabels()
                                 .at(node.GetId().value)
                                 .GetLeafSet()
                                 ->GetClades();
        if (not merge_.ContainsLeafset(clades_difference(clades, src_clades))) {
          ++new_nodes_count;
        }
        curr_node = curr_node->parent;
        if (curr_node == nullptr) {
          break;
        }
      }

      curr_node = move.dst;
      while (not(curr_node->node_id == lca_id.value)) {
        MADAG::NodeView node = merge_.GetResult().Get(NodeId{curr_node->node_id});
        const auto& clades = merge_.GetResultNodeLabels()
                                 .at(node.GetId().value)
                                 .GetLeafSet()
                                 ->GetClades();
        if (not merge_.ContainsLeafset(clades_union(clades, dst_clades))) {
          ++new_nodes_count;
        }
        curr_node = curr_node->parent;
        if (curr_node == nullptr) {
          break;
        }
      }
    }

    move.score_change = move_score_coeffs_.second * move.score_change -
                        move_score_coeffs_.first * new_nodes_count;
    return move.score_change <= 0;
  }

 private:
  const Merge& merge_;
  MADAG sample_;
  const std::vector<NodeId>& sample_dag_ids_;
  const std::pair<int, int> move_score_coeffs_;
};

int main(int argc, char** argv) try {
  Arguments args = GetArguments(argc, argv);
  int ignored{};
  std::string input_dag_path;
  std::string output_dag_path;
  std::string matoptimize_path = "matOptimize";
  std::string logfile_path = "optimization_log";
  std::string refseq_path;
  bool sample_best_tree = false;
  size_t count = 1;
  int move_coeff_nodes = 1;
  int move_coeff_pscore = 1;

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "-i" or name == "--input") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      input_dag_path = *params.begin();
    } else if (name == "-o" or name == "--output") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      output_dag_path = *params.begin();
    } else if (name == "-m" or name == "--matopt") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      matoptimize_path = *params.begin();
    } else if (name == "-c" or name == "--count") {
      if (params.empty()) {
        std::cerr << "Count not specified.\n";
        Fail();
      }
      count = static_cast<size_t>(ParseNumber(*params.begin()));
    } else if (name == "-l" or name == "--logpath") {
      if (params.empty()) {
        std::cerr << "log path name not specified.\n";
        Fail();
      }
      logfile_path = *params.begin();
    } else if (name == "--move-coeff-pscore") {
      if (params.empty()) {
        std::cerr << "parsimony score move coefficient not specified\n";
        Fail();
      }
      move_coeff_pscore = ParseNumber(*params.begin());
    } else if (name == "--move-coeff-nodes") {
      if (params.empty()) {
        std::cerr << "parsimony score move coefficient not specified\n";
        Fail();
      }
      move_coeff_nodes = ParseNumber(*params.begin());
    } else if (name == "--sample-best-tree") {
      sample_best_tree = true;
    } else if (name == "-r" or name == "--MAT-refseq-file") {
      if (params.empty()) {
        std::cerr << "Mutation annotated tree refsequence fasta path not specified.\n";
        Fail();
      }
      refseq_path = *params.begin();
    } else {
      std::cerr << "Unknown argument.\n";
      Fail();
    }
  }

  if (input_dag_path.empty()) {
    std::cerr << "Path to input DAG not specified.\n";
    Fail();
  }

  if (output_dag_path.empty()) {
    std::cerr << "Path to output DAG not specified.\n";
    Fail();
  }
  std::filesystem::create_directory(logfile_path);
  std::string logfile_name = logfile_path + "/logfile.csv";

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ignored);

  std::ofstream logfile;
  logfile.open(logfile_name);
  logfile << "Iteration\tNTrees\tNNodes\tNEdges\tMaxParsimony\tNTreesMaxParsimony\tWors"
             "tParsimony\tSecondsElapsed";

  MADAGStorage input_dag =
      refseq_path.empty()
          ? LoadDAGFromProtobuf(input_dag_path)
          : LoadTreeFromProtobuf(input_dag_path, LoadReferenceSequence(refseq_path));

  input_dag.View().RecomputeCompactGenomes();
  Merge merge{input_dag.View().GetReferenceSequence()};
  merge.AddDAGs({input_dag.View()});
  std::vector<MADAGStorage> optimized_dags;

  auto start_time = std::chrono::high_resolution_clock::now();
  auto time_elapsed = [&start_time]() {
    auto now = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
  };

  auto logger = [&merge, &logfile, &time_elapsed](size_t iteration) {
    SubtreeWeight<BinaryParsimonyScore> parsimonyscorer{merge.GetResult()};
    SubtreeWeight<MaxBinaryParsimonyScore> maxparsimonyscorer{merge.GetResult()};
    merge.ComputeResultEdgeMutations();
    auto minparsimony =
        parsimonyscorer.ComputeWeightBelow(merge.GetResult().GetRoot(), {});
    auto minparsimonytrees =
        parsimonyscorer.MinWeightCount(merge.GetResult().GetRoot(), {});
    auto maxparsimony =
        maxparsimonyscorer.ComputeWeightBelow(merge.GetResult().GetRoot(), {});
    SubtreeWeight<TreeCount> treecount{merge.GetResult()};
    auto ntrees = treecount.ComputeWeightBelow(merge.GetResult().GetRoot(), {});
    std::cout << "Best parsimony score in DAG: " << minparsimony << "\n";
    std::cout << "Worst parsimony score in DAG: " << maxparsimony << "\n";
    std::cout << "Total trees in DAG: " << ntrees << "\n";
    std::cout << "Optimal trees in DAG: " << minparsimonytrees << "\n";
    logfile << '\n'
            << iteration << '\t' << ntrees << '\t' << merge.GetResult().GetNodesCount()
            << '\t' << merge.GetResult().GetEdgesCount() << '\t' << minparsimony << '\t'
            << minparsimonytrees << '\t' << maxparsimony << '\t' << time_elapsed()
            << std::flush;
  };
  logger(0);

  for (size_t i = 0; i < count; ++i) {
    std::cout << "############ Beginning optimize loop " << std::to_string(i)
              << " #######\n";

    merge.ComputeResultEdgeMutations();
    SubtreeWeight<BinaryParsimonyScore> weight{merge.GetResult()};
    auto [sample, dag_ids] =
        sample_best_tree ? weight.MinWeightSampleTree({}) : weight.SampleTree({});
    check_edge_mutations(sample.View());
    Larch_Move_Found_Callback callback{
        merge, sample.View(), dag_ids, {move_coeff_nodes, move_coeff_pscore}};
    /* StoreTreeToProtobuf(sample.View(), "before_optimize_dag.pb"); */
    MADAGStorage result = optimize_dag_direct(sample.View(), callback);
    optimized_dags.push_back(std::move(result));
    merge.AddDAGs({optimized_dags.back().View()});

    if (i % 10 == 0) {  // NOLINT
      StoreDAGToProtobuf(merge.GetResult(), logfile_path + "/intermediate_dag.pb");
    }
    logger(i + 1);
  }

  std::cout << "new node coefficient: " << move_coeff_nodes << "\n";
  std::cout << "parsimony score coefficient: " << move_coeff_pscore << "\n";
  logfile.close();
  StoreDAGToProtobuf(merge.GetResult(), output_dag_path);

  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << "Uncaught exception: " << e.what() << std::endl;
  std::terminate();
} catch (...) {
  std::abort();
}

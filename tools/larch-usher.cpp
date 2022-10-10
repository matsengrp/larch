#include <cstdlib>
#include <iostream>
#include <fstream>
#include <charconv>

#include <unistd.h>
#include <sys/wait.h>

#include "arguments.hpp"
#include "dag_loader.hpp"
#include "mutation_annotated_dag.hpp"
#include "node.hpp"
#include "range/v3/algorithm/sort.hpp"
#include "range/v3/algorithm/unique.hpp"
#include "range/v3/view/zip.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "subtree_weight.hpp"
#include "weight_accumulator.hpp"
#include "tree_count.hpp"
#include "parsimony_score.hpp"
#include "merge.hpp"
#include <mpi.h>

#include "../deps/usher/src/matOptimize/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"

MADAG optimize_dag_direct(const MADAG& dag, Move_Found_Callback& callback);
[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "larch-usher -i,--input file -o,--output file [-m,--matopt file] "
               "[-c,--count number]\n";
  std::cout << "  -i,--input   Path to input DAG\n";
  std::cout << "  -o,--output  Path to output DAG\n";
  std::cout << "  -m,--matopt  Path to matOptimize executable. Default: matOptimize\n";
  std::cout << "  -l,--logfile  Name for logging csv file. Default: logfile.csv\n";
  std::cout << "  -c,--count   Number of iterations. Default: 1\n";

  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

static size_t ParseNumber(std::string_view str) {
  size_t number;
  if (std::from_chars(str.begin(), str.end(), number).ec == std::errc{}) {
    return number;
  }
  throw std::runtime_error("Invalid number");
}

static void CallMatOptimize(std::string matoptimize_path, std::string input,
                            std::string output) {
  pid_t pid = fork();
  if (pid == 0) {
    if (execl(matoptimize_path.c_str(), "matOptimize", "-i", input.c_str(), "-o",
              output.c_str(), "-T", "1", "-n", nullptr) == -1) {
      throw std::runtime_error("Exec failed");
    }
  } else if (pid > 0) {
    int status;
    if (wait(&status) == -1) {
      throw std::runtime_error("Wait failed");
    }
    if (not WIFEXITED(status) or WEXITSTATUS(status) != EXIT_SUCCESS) {
      throw std::runtime_error("Child process failed");
    }
  } else {
    throw std::runtime_error("Fork failed");
  }
}

void check_edge_mutations(const MADAG& madag);

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

struct Larch_Move_Found_Callback : public Move_Found_Callback {
  Larch_Move_Found_Callback(const Merge& merge, const MADAG& sample,
                            const std::vector<NodeId>& sample_dag_ids)
      : merge_{merge}, sample_{sample}, sample_dag_ids_{sample_dag_ids} {}
  bool operator()(Profitable_Moves move, int best_score_change,
                  std::vector<Node_With_Major_Allele_Set_Change>&
                      node_with_major_allele_set_change) override {
    auto& src_clades = merge_.GetResultNodeLabels()
                           .at(sample_dag_ids_.at(move.src->node_id).value)
                           .GetLeafSet()
                           ->GetClades();
    auto& dst_clades = merge_.GetResultNodeLabels()
                           .at(sample_dag_ids_.at(move.dst->node_id).value)
                           .GetLeafSet()
                           ->GetClades();

    return not merge_.ContainsLeafset(clades_union(src_clades, dst_clades)) or
           not merge_.ContainsLeafset(clades_difference(src_clades, dst_clades));
  }
  const Merge& merge_;
  const MADAG& sample_;
  const std::vector<NodeId>& sample_dag_ids_;
};

int main(int argc, char** argv) {
  Arguments args = GetArguments(argc, argv);
  int ignored;
  std::string input_dag_path = "";
  std::string output_dag_path = "";
  std::string matoptimize_path = "matOptimize";
  std::string logfile_path = "logfile.csv";
  size_t count = 1;

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
      count = ParseNumber(*params.begin());
    } else if (name == "-l" or name == "--logfile") {
      if (params.empty()) {
        std::cerr << "Logfile name not specified.\n";
        Fail();
      }
      logfile_path = *params.begin();
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

  auto init_result = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ignored);

  std::ofstream logfile;
  logfile.open(logfile_path);
  logfile << "Iteration\tNTrees\tMaxParsimony\tNTreesMaxParsimony";


  MADAG input_dag = LoadDAGFromProtobuf(input_dag_path);
  Merge merge{input_dag.GetReferenceSequence()};
  merge.AddDAGs({input_dag});
  std::vector<MADAG> optimized_dags;

  auto logger = [&merge, &logfile](size_t iteration) {
      SubtreeWeight<WeightAccumulator<ParsimonyScore>> weightcounter{merge.GetResult()};
      merge.ComputeResultEdgeMutations();
      auto parsimonyscores = weightcounter.ComputeWeightBelow(merge.GetResult().GetDAG().GetRoot(), {});

      std::cout << "Parsimony scores of trees in DAG: "
                << parsimonyscores
                << "\n";

      SubtreeWeight<TreeCount> treecount{merge.GetResult()};
      auto ntrees = treecount.ComputeWeightBelow(merge.GetResult().GetDAG().GetRoot(), {});
      std::cout << "Total trees in DAG: "
                << ntrees
                << "\n";
      logfile << '\n'
              << iteration << '\t'
              << ntrees << '\t'
              << parsimonyscores.GetWeights().begin()->first << '\t'
              << parsimonyscores.GetWeights().begin()->second << '\t';
  };
  logger(0);

  for (size_t i = 0; i < count; ++i) {
    std::cout << "############ Beginning optimize loop " << std::to_string(i)
              << " #######\n";
    merge.ComputeResultEdgeMutations();
    SubtreeWeight<ParsimonyScore> weight{merge.GetResult()};
    auto [sample, dag_ids] = weight.SampleTree({});
    check_edge_mutations(sample);
    MADAG result;
    Larch_Move_Found_Callback callback{merge, sample, dag_ids};
    result = optimize_dag_direct(sample, callback);
    optimized_dags.push_back(std::move(result));
    merge.AddDAGs({optimized_dags.back()});

    logger(i + 1);
  }

  logfile.close();
  StoreDAGToProtobuf(merge.GetResult().GetDAG(),
                     merge.GetResult().GetReferenceSequence(),
                     merge.GetResult().GetEdgeMutations(), output_dag_path);

  return EXIT_SUCCESS;
}

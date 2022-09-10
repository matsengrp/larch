#include <cstdlib>
#include <iostream>
#include <charconv>

#include <unistd.h>
#include <sys/wait.h>

#include "arguments.hpp"
#include "dag_loader.hpp"
#include "subtree_weight.hpp"
#include "parsimony_score.hpp"
#include "merge.hpp"

void assign_sample_ids(MADAG &dag1, MADAG &dag2) {
    std::vector<CompactGenome> &dag1_cgs = dag1.GetCompactGenomes();
    std::vector<CompactGenome> &dag2_cgs = dag2.GetCompactGenomes();

    if (dag1_cgs.size() > 0 and dag2_cgs.size() > 0) {
        for (auto node1 : dag1.GetDAG().GetNodes()) {
            if (node1.IsLeaf()) {
                CompactGenome cg1 = dag1_cgs.at(node1.GetId().value).Copy();
                for (auto node2 : dag2.GetDAG().GetNodes()) {
                     if (node2.IsLeaf()  and 
                         (cg1 == dag2_cgs.at(node2.GetId().value).Copy())) {
                         node2.SetSampleId(node1.GetSampleId());
                     }
                }
            }
        }
    } else {
        throw std::runtime_error("could not assign sample ids when one DAG has no compact genome");
    }
}

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "larch-usher -i,--input file -o,--output file [-m,--matopt file] "
               "[-c,--count number]\n";
  std::cout << "  -i,--input   Path to input DAG\n";
  std::cout << "  -o,--output  Path to output DAG\n";
  std::cout << "  -m,--matopt  Path to matOptimize executable. Default: matOptimize\n";
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

int main(int argc, char** argv) {
  Arguments args = GetArguments(argc, argv);

  std::string input_dag_path = "";
  std::string output_dag_path = "";
  std::string matoptimize_path = "matOptimize";
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

  MADAG input_dag = LoadDAGFromProtobuf(input_dag_path);
  Merge merge{input_dag.GetReferenceSequence()};
  merge.AddDAGs({input_dag});
  std::vector<MADAG> optimized_dags;

  for (size_t i = 0; i < count; ++i) {
    merge.GetResult().GetEdgeMutations() = merge.ComputeResultEdgeMutations();
    SubtreeWeight<ParsimonyScore> weight{merge.GetResult()};
    MADAG sample = weight.SampleTree({});
    input_dag.GetCompactGenomes() = input_dag.ComputeCompactGenomes(input_dag.GetReferenceSequence());
    sample.GetCompactGenomes() = sample.ComputeCompactGenomes(sample.GetReferenceSequence());
    assign_sample_ids(input_dag, sample);

    sample.GetCompactGenomes() = sample.ComputeCompactGenomes(input_dag.GetReferenceSequence());
    std::cout << " got sample tree and computed its compact genome.\n" << std::flush;

    StoreTreeToProtobuf(sample, "sampled_tree.pb");
    std::cout << "wrote sample tree " << std::to_string(i) << "\n" << std::flush;
    CallMatOptimize(matoptimize_path, "sampled_tree.pb", "optimized_tree.pb");

    std::cout << "ran matOptimize on the sample tree, and now we will try to recompute its compact genome\n" << std::flush;
    sample.GetCompactGenomes() = sample.ComputeCompactGenomes(input_dag.GetReferenceSequence());
    std::cout << "...computed its compact genome successfully.\n" << std::flush;

    MADAG sample1 = LoadTreeFromProtobuf("sampled_tree.pb");
    sample1.GetReferenceSequence() = input_dag.GetReferenceSequence();
    sample1.GetCompactGenomes() = sample1.ComputeCompactGenomes(input_dag.GetReferenceSequence());
    std::cout << " got sample tree and computed its compact genome.\n" << std::flush;

    MADAG mo_dag = LoadTreeFromProtobuf("./optimized_tree.pb");
    std::cout << "was able to read the matOptimize tree back from protobuf.\n" << std::flush;
    mo_dag.GetReferenceSequence() = input_dag.GetReferenceSequence();
    mo_dag.GetCompactGenomes() = mo_dag.ComputeCompactGenomes(mo_dag.GetReferenceSequence());
    std::cout << "...computed its compact genome successfully.\n" << std::flush;

    optimized_dags.push_back(LoadTreeFromProtobuf("optimized_tree.pb"));
    merge.AddDAGs({optimized_dags.back()});
  }

  StoreDAGToProtobuf(merge.GetResult().GetDAG(), merge.GetReferenceSequence(),
                     merge.ComputeResultEdgeMutations(), output_dag_path);

  return EXIT_SUCCESS;
}


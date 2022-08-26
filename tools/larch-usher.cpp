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
              output.c_str(), nullptr) == -1) {
      throw std::runtime_error("Exec failed");
    }
  } else if (pid > 0) {
    int status;
    if (wait(&status) == -1) {
      throw std::runtime_error("Wait failed");
    }
    if (not WIFEXITED(status)) {
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

  for (size_t i = 0; i < count; ++i) {
    merge.GetResult().GetEdgeMutations() = merge.ComputeResultEdgeMutations();
    SubtreeWeight<ParsimonyScore> weight{merge.GetResult()};
    MADAG sample = weight.SampleTree({});

    StoreTreeToProtobuf(sample, "sampled_tree.pb");
    CallMatOptimize(matoptimize_path, "sampled_tree.pb", "optimized_tree.pb");
    MADAG optimized = LoadTreeFromProtobuf("optimized_tree.pb");
    merge.AddDAGs({optimized});
  }

  StoreDAGToProtobuf(merge.GetResult().GetDAG(), merge.GetReferenceSequence(),
                     merge.ComputeResultEdgeMutations(), output_dag_path);

  return EXIT_SUCCESS;
}

#include <iostream>
#include <cstdlib>
#include <random>

#include "arguments.hpp"
#include "usher_optimize.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "optimize -i,--input file\n";
  std::cout << "  -i,--input     Input tree file\n";
  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

static int Optimize(std::string_view input_path) {
  std::cout << "Optimizing: " << input_path << "\n";

  Mutation_Annotated_Tree::Tree tree;
  Mutation_Annotated_Tree::load_mutation_annotated_tree(std::string{input_path}, tree);

  UsherOptimize(tree);

  return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
  InitUsherMPI(argc, argv);

  Arguments args = GetArguments(argc, argv);

  std::string input_path;

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "-i" or name == "--input") {
      if (params.empty()) {
        std::cerr << "Specify input file name.\n";
        Fail();
      }
      input_path = *params.begin();
    }
  }

  if (input_path.empty()) {
    std::cerr << "Specify at input file name.\n";
    Fail();
  }

  return Optimize(input_path);
}

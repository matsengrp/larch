#include <iostream>
#include <cstdlib>
#include <random>

#include "arguments.hpp"
#include "usher_optimize.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "optimize -i,--input file\n";
  std::cout << "  -i,--input     Input file\n";
  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

static int Optimize(std::string_view input_path) {
  std::cout << "Optimizing: " << input_path << "\n";

  Mutation_Annotated_Tree::Tree tree =
      MAT::load_mutation_annotated_tree(std::string{input_path});
  tree.uncondense_leaves();
  tree.populate_ignored_range();

  UsherOptimize(tree);

  return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
  int ignored;
  auto init_result = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ignored);
  if (init_result != MPI_SUCCESS) {
    fprintf(stderr, "MPI init failed\n");
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &process_count);
  fprintf(stderr, "Running with %d processes\n", process_count);

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

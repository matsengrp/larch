#include <cstdlib>
#include <iostream>

#include "arguments.hpp"
#include "merge.hpp"
#include "history_dag_loader.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "dag_diff -p,--proto file -j,--json file"
               "[-o,--output filename]\n";
  std::cout << "  -p,--proto     Input protobuf DAG filename\n";
  std::cout << "  -j,--json      Input json DAG filename\n";

  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

static int TakeDiff(std::string_view proto_filename, std::string_view json_filename) {
  std::string lhs_reference_sequence;
  std::vector<std::vector<Mutations>> lhs_mutations;
  std::vector<HistoryDAG> lhs_trees;
  
  lhs_mutations.push_back({});
  lhs_trees.push_back(LoadHistoryDAGFromProtobufGZ(proto_filename, lhs_reference_sequence, lhs_mutations.at(0)));
  Merge lhs_merge{lhs_reference_sequence};
  std::vector<std::reference_wrapper<const HistoryDAG>> lhs_tree_refs{lhs_trees.begin(),
                                                                  lhs_trees.end()};
  lhs_merge.AddTrees(lhs_tree_refs, lhs_mutations);

  return EXIT_SUCCESS;
}

int main(int argc, char** argv) {
  Arguments args = GetArguments(argc, argv);

  std::string_view proto_filename;
  std::string_view json_filename;

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "-p" or name == "--proto") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      proto_filename = *params.begin();
    } else if (name == "-j" or name == "--json") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      json_filename = *params.begin();
    } else {
      std::cerr << "Unknown argument.\n";
      Fail();
    }
  }

  if (proto_filename.empty() or json_filename.empty()) {
    std::cerr << "Specify two input file names.\n";
    Fail();
  }

  return TakeDiff(proto_filename, json_filename);
}
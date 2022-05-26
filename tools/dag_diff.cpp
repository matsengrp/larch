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

// static Merge LoadProtobuf(std::string_view path) {
//   std::vector<Mutations> mutations;
//   std::string reference_sequence;

//   HistoryDAG dag = LoadHistoryDAGFromJsonGZ(path, reference_sequence, mutations);
// }

// static Merge LoadJson(std::string_view path) {

// }

static int TakeDiff(std::string_view proto_filename, std::string_view json_filename) {
  //   Merge lhs = LoadProtobuf(proto_filename);
  //   Merge rhs = LoadJson(json_filename);
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

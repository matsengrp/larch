#include <cstdlib>
#include <iostream>

#include "arguments.hpp"
#include "larch/dag_loader.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "dag2dot -t,--tree-pb file\n";
  std::cout << "dag2dot [-c,--cgs] -d,--dag-pb file\n";
  std::cout << "dag2dot -j,--dag-json file\n";
  std::cout << "  -t,--tree-pb   Input protobuf tree filename\n";
  std::cout << "  -d,--dag-pb    Input protobuf DAG filename\n";
  std::cout << "  -j,--dag-json  Input json DAG filename\n";

  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

enum class InputType { TreePB, DagPB, DagJson };

static MADAGStorage Load(InputType type, std::string_view path) {
  switch (type) {
    case InputType::TreePB:
      return LoadTreeFromProtobuf(path, "");
    case InputType::DagPB: {
      MADAGStorage dag = LoadDAGFromProtobuf(path);
      dag.View().RecomputeCompactGenomes();
      return dag;
    }
    case InputType::DagJson:
      return LoadDAGFromJson(path);
    default:
      std::cerr << "Unknown input type\n";
      Fail();
  }
}

int main(int argc, char** argv) {
  Arguments args = GetArguments(argc, argv);

  InputType type = InputType::TreePB;
  std::string path;

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "-t" or name == "--tree-pb") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      type = InputType::TreePB;
      path = *params.begin();
    } else if (name == "-d" or name == "--dag-pb") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      type = InputType::DagPB;
      path = *params.begin();
    } else if (name == "-j" or name == "--dag-json") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      type = InputType::DagJson;
      path = *params.begin();
    } else {
      std::cerr << "Unknown argument.\n";
      Fail();
    }
  }

  if (path.empty()) {
    std::cerr << "Filename not specified.\n";
    Fail();
  }

  MADAGStorage dag = Load(type, path);

  MADAGToDOT(dag.View(), std::cout);

  return EXIT_SUCCESS;
}

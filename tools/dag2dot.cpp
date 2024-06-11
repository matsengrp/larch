#include <cstdlib>
#include <iostream>

#include "arguments.hpp"
#include "larch/dag_loader.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "dag2dot -i,--input file\n";
  std::cout << "  -i,--input     Input filename (format inferred)\n";
  std::cout << "  -t,--tree-pb   Input protobuf tree filename\n";
  std::cout << "  -d,--dag-pb    Input protobuf DAG filename\n";
  std::cout << "  -j,--dag-json  Input json DAG filename\n";
  std::cout << "  -b,--dagbin    Input dagbin DAG filename\n";
  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";
  std::exit(EXIT_FAILURE);
}

int main(int argc, char** argv) try {
  Arguments args = GetArguments(argc, argv);

  std::string input_path;
  FileFormat input_format = FileFormat::Infer;

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "-i" or name == "--input") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      input_path = *params.begin();
    } else if (name == "-t" or name == "--tree-pb") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      input_format = FileFormat::ProtobufTree;
      input_path = *params.begin();
    } else if (name == "-d" or name == "--dag-pb") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      input_format = FileFormat::ProtobufDAG;
      input_path = *params.begin();
    } else if (name == "-j" or name == "--dag-json") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      input_format = FileFormat::JsonDAG;
      input_path = *params.begin();
    } else if (name == "-b" or name == "--dagbin") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      input_format = FileFormat::JsonDAG;
      input_path = *params.begin();
    } else {
      std::cerr << "Unknown argument '" << name << "'.\n";
      Fail();
    }
  }

  if (input_path.empty()) {
    std::cerr << "Input path not specified.\n";
    Fail();
  }
  if (input_format == FileFormat::Infer) {
    input_format = InferFileFormat(input_path);
  }

  auto dag = LoadDAG(input_path, input_format);
  MADAGToDOT(dag.View(), std::cout);

  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << "Uncaught exception: " << e.what() << std::endl;
  std::terminate();
} catch (...) {
  std::abort();
}

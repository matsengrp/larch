#include <cstdlib>
#include <iostream>

#include "arguments.hpp"
#include "larch/dag_loader.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "dag2dot -i,--input file\n";
  std::cout << "  -i,--input       Path to input Tree/DAG\n";
  std::cout << "  -o,--output      Path to output DOT file (default: \n";
  std::cout << "  --input-format   Input file format (default: inferred)\n";
  std::cout << "  --dag/--tree     Specify whether input is a DAG or Tree\n";
  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";
  std::exit(EXIT_FAILURE);
}

int main(int argc, char** argv) try {
  Arguments args = GetArguments(argc, argv);

  std::string input_path;
  std::string output_path;
  FileFormat input_format = FileFormat::Infer;
  std::optional<bool> is_input_dag = std::nullopt;

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "-i" or name == "--input") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      input_path = *params.begin();
    } else if (name == "-o" or name == "--output") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      output_path = *params.begin();
    } else if (name == "--input-format") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      input_format = FileFormat::ProtobufTree;
      input_path = *params.begin();
    } else if (name == "--dag" or name == "--tree") {
      is_input_dag = (name == "--dag");
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
  if (input_format == FileFormat::Protobuf) {
    if (is_input_dag.has_value()) {
      input_format =
          is_input_dag.value() ? FileFormat::ProtobufDAG : FileFormat::ProtobufTree;
    } else {
      std::cerr
          << "ERROR: When using *.pb protobuf files, please specify --dag or --tree."
          << std::endl;
      Fail();
    }
  }

  auto dag = LoadDAG(input_path, input_format);

  if (output_path.empty()) {
    MADAGToDOT(dag.View(), std::cout);
  } else {
    std::ofstream outfile(output_path);
    MADAGToDOT(dag.View(), outfile);
    outfile.close();
  }

  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << "Uncaught exception: " << e.what() << std::endl;
  std::terminate();
} catch (...) {
  std::abort();
}

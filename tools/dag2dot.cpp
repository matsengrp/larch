#include <cstdlib>
#include <iostream>

#include "tools_common.hpp"
#include "larch/dag_loader.hpp"

[[noreturn]] static void Usage() {
  std::string program_desc = "dag2dot: tool for converting DAG/tree file to DOT file";

  std::vector<std::string> usage_examples = {
      {"dag2dot -i,--input FILE [-o,--output FILE]"}};

  std::vector<std::pair<std::string, std::string>> flag_desc_pairs = {
      {"-i,--input FILE", "Path to input DAG/Tree file (REQUIRED)"},
      {"-o,--output FILE", "Path to output DOT file (default: DOT written to stdout)"},
      {"--input-format ENUM",
       "Specify input file format (default: inferred) \n"
       "[dagbin, dag-pb, tree-pb, dag-json]"},
      {"--dag/--tree", "Specify whether protobuf input is a DAG or Tree"}};

  std::cout << FormatUsage(program_desc, usage_examples, flag_desc_pairs);

  std::exit(EXIT_SUCCESS);
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
      ParseOption(name, params, input_path, 1);
    } else if (name == "-o" or name == "--output") {
      ParseOption(name, params, output_path, 1);
    } else if (name == "--input-format") {
      std::string temp;
      ParseOption(name, params, temp, 1);
      input_format = InferFileFormat(temp);
    } else if (name == "--dag" or name == "--tree") {
      ParseOption<false>(name, params, is_input_dag, 0);
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

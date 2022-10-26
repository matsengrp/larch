#include <cstdlib>
#include <iostream>
#include <fstream>

#include "arguments.hpp"
#include "dag_loader.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "mat2dag -t,--tree-pb input_file -r --refseq refseq_file -o --output-pb output_file\n";
  std::cout << "  -t,--tree-pb              Input protobuf tree file\n";
  std::cout << "  -r,--refseq               Reference Sequence file\n";
  std::cout << "  -o,--output-pb [OPTIONAL] Output protobuf dag file\n";

  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

bool check_valid(MADAG &madag) {
  // check edge mutations have correct size
  if (madag.GetDAG().GetEdgesCount() != madag.GetEdgeMutations().size()) {
    std::cout << "failed with edge mutations size\n" << std::flush;
    return false;
  }

  // check dag's graph structure is tree shaped (since it comes from a MAT pb)
  if (madag.GetDAG().GetNodesCount() != madag.GetDAG().GetEdgesCount() + 1) {
    std::cout << "failed with tree shape\n" << std::flush;
    return false;
  }

  // check node labels, as given by compact genome, are unique
  madag.RemoveCompactGenomes();
  madag.RecomputeCompactGenomes();

  //madag.AssertUA();
  return true;
}

int main(int argc, char** argv) {
  Arguments args = GetArguments(argc, argv);

  std::string pb_path = "";
  std::string seq_path = "";
  std::string output_path = "";

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "-t" or name == "--tree-pb") {
      if (params.empty()) {
        std::cerr << "MAT protobuf filename not specified.\n";
        Fail();
      }
      pb_path = *params.begin();
    } else if (name == "-r" or name == "--refseq") {
      if (params.empty()) {
        std::cerr << "Reference sequence filename not specified.\n";
        Fail();
      }
      seq_path = *params.begin();
    } else if (name == "-o" or name == "--output-pb") {
      if (params.empty()) {
        std::cerr << "Output filename not specified.\n";
        Fail();
      }
      output_path = *params.begin();
    } else {
      std::cerr << "Unknown argument.\n";
      Fail();
    }
  }

  if (pb_path.empty()) {
    std::cerr << "MAT protobuf filename not specified.\n";
    Fail();
  }
  if (seq_path.empty()) {
    std::cerr << "Reference sequence filename not specified.\n";
    Fail();
  }

  MADAG dag;
  std::string reference_sequence = "";
  std::fstream file;
  file.open(seq_path);
  while (file >> reference_sequence) {}
  dag = LoadTreeFromProtobuf(pb_path, reference_sequence);

  Assert(check_valid(dag));

  if (!output_path.empty()) {
    StoreDAGToProtobuf(dag.GetDAG(), dag.GetReferenceSequence(), dag.GetEdgeMutations(), output_path);
  }
  return EXIT_SUCCESS;
}

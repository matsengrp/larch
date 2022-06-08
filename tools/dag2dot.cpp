#include <cstdlib>
#include <iostream>

#include "arguments.hpp"
#include "dag_loader.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "dag2dot -t,--tree-pb file\n";
  std::cout << "dag2dot [-c,--cgs] -d,--dag-pb file\n";
  std::cout << "dag2dot [-m,--muts] -j,--dag-json file\n";
  std::cout << "  -t,--tree-pb   Input protobuf tree filename\n";
  std::cout << "  -d,--dag-pb    Input protobuf DAG filename\n";
  std::cout << "  -j,--dag-json  Input json DAG filename\n";
  std::cout << "  -c,--cgs       Compute compact genomes\n";
  std::cout << "  -m,--muts      Compute edge mutations\n";

  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

enum class InputType { TreePB, DagPB, DagJson };

int main(int argc, char** argv) {
  Arguments args = GetArguments(argc, argv);

  InputType type = InputType::TreePB;
  std::string path = "";
  bool cgs = false;
  bool muts = false;

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
    } else if (name == "-c" or name == "--cgs") {
      cgs = true;
    } else if (name == "-m" or name == "--muts") {
      muts = true;
    } else {
      std::cerr << "Unknown argument.\n";
      Fail();
    }
  }

  if (path.empty()) {
    std::cerr << "Filename not specified.\n";
    Fail();
  }

  MADAG dag;
  switch (type) {
    case InputType::TreePB:
      dag = LoadTreeFromProtobuf(path);
      dag.dag.ReindexPreOrder();
      EdgeMutationsDAGToDOT(dag, std::cout);
      break;
    case InputType::DagPB:
      dag = LoadDAGFromProtobuf(path);
      dag.dag.ReindexPreOrder();
      if (cgs) {
        Assert(not dag.reference_sequence.empty());
        CompactGenomeDAGToDOT(dag, dag.ComputeCompactGenomesDAG(dag.reference_sequence),
                              std::cout);
      } else {
        EdgeMutationsDAGToDOT(dag, std::cout);
      }
      break;
    case InputType::DagJson:
      dag = LoadDAGFromJson(path);
      std::map<NodeId, NodeId> index = dag.dag.ReindexPreOrder();
      std::vector<CompactGenome> compact_genomes = LoadCompactGenomesJson(path);
      {
        std::vector<CompactGenome> reindexed_cgs;
        reindexed_cgs.reserve(compact_genomes.size());
        for (auto [prev, curr] : index) {
          reindexed_cgs.emplace_back(std::move(compact_genomes.at(curr.value)));
        }
        compact_genomes = std::move(reindexed_cgs);
      }
      if (muts) {
        dag.edge_mutations.resize(dag.dag.GetEdges().size());
        for (Edge edge : dag.dag.GetEdges()) {
          auto [parent, child] = edge.GetNodeIds();
          dag.edge_mutations.at(edge.GetId().value) = CompactGenome::ToEdgeMutations(
              dag.reference_sequence, compact_genomes.at(parent.value),
              compact_genomes.at(child.value));
        }
        EdgeMutationsDAGToDOT(dag, std::cout);
      } else {
        CompactGenomeDAGToDOT(dag, compact_genomes, std::cout);
      }
      break;
  }

  return EXIT_SUCCESS;
}

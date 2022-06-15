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

MADAG MakeSyntheticDAG() {
  MADAG result;

  result.GetReferenceSequence() = "ACGTACGT";

  result.GetDAG().InitializeNodes(31);
  result.GetDAG().AppendEdge({0}, {1}, {0});    // 0
  result.GetDAG().AppendEdge({1}, {2}, {0});    // 1
  result.GetDAG().AppendEdge({1}, {3}, {1});    // 2
  result.GetDAG().AppendEdge({2}, {4}, {0});    // 3
  result.GetDAG().AppendEdge({2}, {5}, {1});    // 4
  result.GetDAG().AppendEdge({2}, {6}, {2});    // 5
  result.GetDAG().AppendEdge({4}, {7}, {0});    // 6
  result.GetDAG().AppendEdge({4}, {8}, {1});    // 7
  result.GetDAG().AppendEdge({5}, {16}, {0});   // 8
  result.GetDAG().AppendEdge({5}, {9}, {1});    // 9
  result.GetDAG().AppendEdge({5}, {11}, {2});   // 10
  result.GetDAG().AppendEdge({6}, {10}, {0});   // 11
  result.GetDAG().AppendEdge({6}, {11}, {1});   // 12
  result.GetDAG().AppendEdge({7}, {12}, {0});   // 13
  result.GetDAG().AppendEdge({7}, {13}, {1});   // 14
  result.GetDAG().AppendEdge({8}, {14}, {0});   // 15
  result.GetDAG().AppendEdge({8}, {15}, {1});   // 16
  result.GetDAG().AppendEdge({8}, {16}, {2});   // 17
  result.GetDAG().AppendEdge({9}, {17}, {0});   // 18
  result.GetDAG().AppendEdge({9}, {18}, {1});   // 19
  result.GetDAG().AppendEdge({9}, {19}, {2});   // 20
  result.GetDAG().AppendEdge({11}, {20}, {0});  // 21
  result.GetDAG().AppendEdge({11}, {21}, {1});  // 22
  result.GetDAG().AppendEdge({15}, {22}, {0});  // 23
  result.GetDAG().AppendEdge({15}, {23}, {1});  // 24
  result.GetDAG().AppendEdge({18}, {24}, {0});  // 25
  result.GetDAG().AppendEdge({18}, {25}, {1});  // 26
  result.GetDAG().AppendEdge({20}, {25}, {0});  // 27
  result.GetDAG().AppendEdge({20}, {26}, {1});  // 28
  result.GetDAG().AppendEdge({22}, {27}, {0});  // 29
  result.GetDAG().AppendEdge({22}, {28}, {1});  // 30
  result.GetDAG().AppendEdge({28}, {29}, {0});  // 31
  result.GetDAG().AppendEdge({28}, {30}, {1});  // 32

  result.GetDAG().BuildConnections();

  result.GetEdgeMutations().resize(result.GetDAG().GetEdgesCount());

  result.GetEdgeMutations().at(1)[{1}] = {'A', 'G'};
  result.GetEdgeMutations().at(2)[{3}] = {'G', 'C'};
  result.GetEdgeMutations().at(3)[{2}] = {'C', 'T'};
  result.GetEdgeMutations().at(4)[{2}] = {'C', 'A'};
  result.GetEdgeMutations().at(5)[{2}] = {'C', 'G'};
  result.GetEdgeMutations().at(5)[{4}] = {'T', 'A'};
  result.GetEdgeMutations().at(6)[{5}] = {'A', 'G'};
  result.GetEdgeMutations().at(8)[{2}] = {'A', 'G'};
  result.GetEdgeMutations().at(8)[{6}] = {'C', 'A'};
  result.GetEdgeMutations().at(9)[{8}] = {'T', 'C'};
  result.GetEdgeMutations().at(10)[{2}] = {'A', 'C'};
  result.GetEdgeMutations().at(10)[{3}] = {'G', 'A'};
  result.GetEdgeMutations().at(11)[{4}] = {'A', 'T'};
  result.GetEdgeMutations().at(12)[{2}] = {'G', 'C'};
  result.GetEdgeMutations().at(13)[{1}] = {'G', 'C'};
  result.GetEdgeMutations().at(13)[{2}] = {'T', 'A'};
  result.GetEdgeMutations().at(14)[{4}] = {'T', 'C'};
  result.GetEdgeMutations().at(15)[{3}] = {'G', 'A'};
  result.GetEdgeMutations().at(16)[{6}] = {'C', 'T'};
  result.GetEdgeMutations().at(17)[{2}] = {'T', 'G'};
  result.GetEdgeMutations().at(17)[{4}] = {'T', 'C'};
  result.GetEdgeMutations().at(18)[{8}] = {'C', 'G'};
  result.GetEdgeMutations().at(19)[{8}] = {'C', 'A'};
  result.GetEdgeMutations().at(20)[{1}] = {'G', 'T'};
  result.GetEdgeMutations().at(23)[{2}] = {'T', 'C'};
  result.GetEdgeMutations().at(25)[{1}] = {'G', 'C'};
  result.GetEdgeMutations().at(26)[{5}] = {'A', 'C'};
  result.GetEdgeMutations().at(26)[{8}] = {'A', 'G'};
  result.GetEdgeMutations().at(27)[{5}] = {'A', 'C'};
  result.GetEdgeMutations().at(28)[{3}] = {'A', 'T'};
  result.GetEdgeMutations().at(30)[{3}] = {'G', 'C'};
  result.GetEdgeMutations().at(31)[{6}] = {'T', 'A'};

  result.GetCompactGenomes() =
      result.ComputeCompactGenomes(result.GetReferenceSequence());

  return result;
}

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
      break;
    case InputType::DagPB:
      dag = LoadDAGFromProtobuf(path);
      if (cgs) {
        Assert(not dag.GetReferenceSequence().empty());
        dag.GetCompactGenomes() = dag.ComputeCompactGenomes(dag.GetReferenceSequence());
      }
      break;
    case InputType::DagJson:
      dag = LoadDAGFromJson(path);
      if (muts) {
        dag.GetEdgeMutations().resize(dag.GetDAG().GetEdgesCount());
        for (Edge edge : dag.GetDAG().GetEdges()) {
          auto [parent, child] = edge.GetNodeIds();
          dag.GetEdgeMutations().at(edge.GetId().value) =
              CompactGenome::ToEdgeMutations(dag.GetReferenceSequence(),
                                             dag.GetCompactGenomes().at(parent.value),
                                             dag.GetCompactGenomes().at(child.value));
        }
      }
      break;
  }
  MADAGToDOT(dag, std::cout);

  return EXIT_SUCCESS;
}

#include "synthetic_dags.hpp"

MADAG MakeSyntheticDAG() {
  MADAG result;

  result.reference_sequence = "ACGTACGT";

  result.dag.InitializeNodes(31);
  result.dag.AppendEdge({0}, {1}, {0});    // 0
  result.dag.AppendEdge({1}, {2}, {0});    // 1
  result.dag.AppendEdge({1}, {3}, {1});    // 2
  result.dag.AppendEdge({2}, {4}, {0});    // 3
  result.dag.AppendEdge({2}, {5}, {1});    // 4
  result.dag.AppendEdge({2}, {6}, {2});    // 5
  result.dag.AppendEdge({4}, {7}, {0});    // 6
  result.dag.AppendEdge({4}, {8}, {1});    // 7
  result.dag.AppendEdge({5}, {16}, {0});   // 8
  result.dag.AppendEdge({5}, {9}, {1});    // 9
  result.dag.AppendEdge({5}, {11}, {2});   // 10
  result.dag.AppendEdge({6}, {10}, {0});   // 11
  result.dag.AppendEdge({6}, {11}, {1});   // 12
  result.dag.AppendEdge({7}, {12}, {0});   // 13
  result.dag.AppendEdge({7}, {13}, {1});   // 14
  result.dag.AppendEdge({8}, {14}, {0});   // 15
  result.dag.AppendEdge({8}, {15}, {1});   // 16
  result.dag.AppendEdge({8}, {16}, {2});   // 17
  result.dag.AppendEdge({9}, {17}, {0});   // 18
  result.dag.AppendEdge({9}, {18}, {1});   // 19
  result.dag.AppendEdge({9}, {19}, {2});   // 20
  result.dag.AppendEdge({11}, {20}, {0});  // 21
  result.dag.AppendEdge({11}, {21}, {1});  // 22
  result.dag.AppendEdge({15}, {22}, {0});  // 23
  result.dag.AppendEdge({15}, {23}, {1});  // 24
  result.dag.AppendEdge({18}, {24}, {0});  // 25
  result.dag.AppendEdge({18}, {25}, {1});  // 26
  result.dag.AppendEdge({20}, {25}, {0});  // 27
  result.dag.AppendEdge({20}, {26}, {1});  // 28
  result.dag.AppendEdge({22}, {27}, {0});  // 29
  result.dag.AppendEdge({22}, {28}, {1});  // 30
  result.dag.AppendEdge({28}, {29}, {0});  // 31
  result.dag.AppendEdge({28}, {30}, {1});  // 32

  result.dag.BuildConnections();

  result.edge_mutations.resize(result.dag.GetEdgesCount());

  result.edge_mutations.at(1)[{1}] = {'A', 'G'};
  result.edge_mutations.at(2)[{3}] = {'G', 'C'};
  result.edge_mutations.at(3)[{2}] = {'C', 'T'};
  result.edge_mutations.at(4)[{2}] = {'C', 'A'};
  result.edge_mutations.at(5)[{2}] = {'C', 'G'};
  result.edge_mutations.at(5)[{4}] = {'T', 'A'};
  result.edge_mutations.at(6)[{5}] = {'A', 'G'};
  result.edge_mutations.at(8)[{2}] = {'A', 'G'};
  result.edge_mutations.at(8)[{6}] = {'C', 'A'};
  result.edge_mutations.at(9)[{8}] = {'T', 'C'};
  result.edge_mutations.at(10)[{2}] = {'A', 'C'};
  result.edge_mutations.at(10)[{3}] = {'G', 'A'};
  result.edge_mutations.at(11)[{4}] = {'A', 'T'};
  result.edge_mutations.at(12)[{2}] = {'G', 'C'};
  result.edge_mutations.at(13)[{1}] = {'G', 'C'};
  result.edge_mutations.at(13)[{2}] = {'T', 'A'};
  result.edge_mutations.at(14)[{4}] = {'T', 'C'};
  result.edge_mutations.at(15)[{3}] = {'G', 'A'};
  result.edge_mutations.at(16)[{6}] = {'C', 'T'};
  result.edge_mutations.at(17)[{2}] = {'T', 'G'};
  result.edge_mutations.at(17)[{4}] = {'T', 'C'};
  result.edge_mutations.at(18)[{8}] = {'C', 'G'};
  result.edge_mutations.at(19)[{8}] = {'C', 'A'};
  result.edge_mutations.at(20)[{1}] = {'G', 'T'};
  result.edge_mutations.at(23)[{2}] = {'T', 'C'};
  result.edge_mutations.at(25)[{1}] = {'G', 'C'};
  result.edge_mutations.at(26)[{5}] = {'A', 'C'};
  result.edge_mutations.at(26)[{8}] = {'A', 'G'};
  result.edge_mutations.at(27)[{5}] = {'A', 'C'};
  result.edge_mutations.at(28)[{3}] = {'A', 'T'};
  result.edge_mutations.at(30)[{3}] = {'G', 'C'};
  result.edge_mutations.at(31)[{6}] = {'T', 'A'};

  result.compact_genomes = result.ComputeCompactGenomes(result.reference_sequence);

  return result;
}

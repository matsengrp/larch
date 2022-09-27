#include "synthetic_dags.hpp"

MADAG MakeSyntheticDAG() {
  MADAG result;

  result.SetReferenceSequence("ACGTACGT");

  result.InitializeNodes(31);
  result.AppendEdge({0}, {1}, {0});    // 0
  result.AppendEdge({1}, {2}, {0});    // 1
  result.AppendEdge({1}, {3}, {1});    // 2
  result.AppendEdge({2}, {4}, {0});    // 3
  result.AppendEdge({2}, {5}, {1});    // 4
  result.AppendEdge({2}, {6}, {2});    // 5
  result.AppendEdge({4}, {7}, {0});    // 6
  result.AppendEdge({4}, {8}, {1});    // 7
  result.AppendEdge({5}, {16}, {0});   // 8
  result.AppendEdge({5}, {9}, {1});    // 9
  result.AppendEdge({5}, {11}, {2});   // 10
  result.AppendEdge({6}, {10}, {0});   // 11
  result.AppendEdge({6}, {11}, {1});   // 12
  result.AppendEdge({7}, {12}, {0});   // 13
  result.AppendEdge({7}, {13}, {1});   // 14
  result.AppendEdge({8}, {14}, {0});   // 15
  result.AppendEdge({8}, {15}, {1});   // 16
  result.AppendEdge({8}, {16}, {2});   // 17
  result.AppendEdge({9}, {17}, {0});   // 18
  result.AppendEdge({9}, {18}, {1});   // 19
  result.AppendEdge({9}, {19}, {2});   // 20
  result.AppendEdge({11}, {20}, {0});  // 21
  result.AppendEdge({11}, {21}, {1});  // 22
  result.AppendEdge({15}, {22}, {0});  // 23
  result.AppendEdge({15}, {23}, {1});  // 24
  result.AppendEdge({18}, {24}, {0});  // 25
  result.AppendEdge({18}, {25}, {1});  // 26
  result.AppendEdge({20}, {25}, {0});  // 27
  result.AppendEdge({20}, {26}, {1});  // 28
  result.AppendEdge({22}, {27}, {0});  // 29
  result.AppendEdge({22}, {28}, {1});  // 30
  result.AppendEdge({28}, {29}, {0});  // 31
  result.AppendEdge({28}, {30}, {1});  // 32

  result.BuildConnections();

  std::vector<EdgeMutations> edge_mutations;
  edge_mutations.resize(result.GetDAG().GetEdgesCount());

  edge_mutations.at(1)[{1}] = {'A', 'G'};
  edge_mutations.at(2)[{3}] = {'G', 'C'};
  edge_mutations.at(3)[{2}] = {'C', 'T'};
  edge_mutations.at(4)[{2}] = {'C', 'A'};
  edge_mutations.at(5)[{2}] = {'C', 'G'};
  edge_mutations.at(5)[{4}] = {'T', 'A'};
  edge_mutations.at(6)[{5}] = {'A', 'G'};
  edge_mutations.at(8)[{2}] = {'A', 'G'};
  edge_mutations.at(8)[{6}] = {'C', 'A'};
  edge_mutations.at(9)[{8}] = {'T', 'C'};
  edge_mutations.at(10)[{2}] = {'A', 'C'};
  edge_mutations.at(10)[{3}] = {'G', 'A'};
  edge_mutations.at(11)[{4}] = {'A', 'T'};
  edge_mutations.at(12)[{2}] = {'G', 'C'};
  edge_mutations.at(13)[{1}] = {'G', 'C'};
  edge_mutations.at(13)[{2}] = {'T', 'A'};
  edge_mutations.at(14)[{4}] = {'T', 'C'};
  edge_mutations.at(15)[{3}] = {'G', 'A'};
  edge_mutations.at(16)[{6}] = {'C', 'T'};
  edge_mutations.at(17)[{2}] = {'T', 'G'};
  edge_mutations.at(17)[{4}] = {'T', 'C'};
  edge_mutations.at(18)[{8}] = {'C', 'G'};
  edge_mutations.at(19)[{8}] = {'C', 'A'};
  edge_mutations.at(20)[{1}] = {'G', 'T'};
  edge_mutations.at(23)[{2}] = {'T', 'C'};
  edge_mutations.at(25)[{1}] = {'G', 'C'};
  edge_mutations.at(26)[{5}] = {'A', 'C'};
  edge_mutations.at(26)[{8}] = {'A', 'G'};
  edge_mutations.at(27)[{5}] = {'A', 'C'};
  edge_mutations.at(28)[{3}] = {'A', 'T'};
  edge_mutations.at(30)[{3}] = {'G', 'C'};
  edge_mutations.at(31)[{6}] = {'T', 'A'};

  result.SetEdgeMutations(std::move(edge_mutations));

  result.RecomputeCompactGenomes();

  return result;
}

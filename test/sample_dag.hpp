#pragma once

#include "larch/madag/mutation_annotated_dag.hpp"

static auto MakeSampleDAG() {
  MADAGStorage input_storage;
  auto dag = input_storage.View();

  dag.SetReferenceSequence("GAA");

  dag.InitializeNodes(11);

  dag.AddEdge({0}, {0}, {10}, {0});
  dag.AddEdge({1}, {7}, {1}, {0}).GetMutableEdgeMutations()[{1}] = {'T', 'A'};
  dag.AddEdge({2}, {7}, {2}, {1}).GetMutableEdgeMutations()[{1}] = {'T', 'G'};
  dag.AddEdge({3}, {8}, {3}, {0}).GetMutableEdgeMutations()[{1}] = {'C', 'A'};
  dag.AddEdge({4}, {8}, {4}, {1}).GetMutableEdgeMutations()[{1}] = {'C', 'A'};
  dag.AddEdge({5}, {9}, {5}, {0}).GetMutableEdgeMutations()[{1}] = {'A', 'C'};
  dag.AddEdge({6}, {9}, {6}, {1}).GetMutableEdgeMutations()[{1}] = {'A', 'T'};
  dag.AddEdge({7}, {8}, {7}, {2}).GetMutableEdgeMutations()[{1}] = {'C', 'T'};
  dag.AddEdge({8}, {10}, {8}, {0}).GetMutableEdgeMutations()[{1}] = {'G', 'C'};
  dag.AddEdge({9}, {10}, {9}, {1}).GetMutableEdgeMutations()[{1}] = {'G', 'A'};

  dag.BuildConnections();

  dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{2}] = {'G', 'C'};
  dag.Get(EdgeId{2}).GetMutableEdgeMutations()[{2}] = {'G', 'T'};
  dag.Get(EdgeId{3}).GetMutableEdgeMutations()[{2}] = {'T', 'G'};
  dag.Get(EdgeId{4}).GetMutableEdgeMutations()[{2}] = {'T', 'C'};
  dag.Get(EdgeId{5}).GetMutableEdgeMutations()[{2}] = {'G', 'T'};
  dag.Get(EdgeId{6}).GetMutableEdgeMutations()[{2}] = {'G', 'C'};
  dag.Get(EdgeId{7}).GetMutableEdgeMutations()[{2}] = {'T', 'G'};
  dag.Get(EdgeId{8}).GetMutableEdgeMutations()[{2}] = {'A', 'T'};
  dag.Get(EdgeId{9}).GetMutableEdgeMutations()[{2}] = {'A', 'G'};

  dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{3}] = {'G', 'C'};
  dag.Get(EdgeId{2}).GetMutableEdgeMutations()[{3}] = {'G', 'T'};
  dag.Get(EdgeId{3}).GetMutableEdgeMutations()[{3}] = {'T', 'G'};
  dag.Get(EdgeId{4}).GetMutableEdgeMutations()[{3}] = {'T', 'G'};
  dag.Get(EdgeId{5}).GetMutableEdgeMutations()[{3}] = {'G', 'T'};
  dag.Get(EdgeId{6}).GetMutableEdgeMutations()[{3}] = {'G', 'C'};
  dag.Get(EdgeId{7}).GetMutableEdgeMutations()[{3}] = {'T', 'G'};
  dag.Get(EdgeId{8}).GetMutableEdgeMutations()[{3}] = {'A', 'C'};
  dag.Get(EdgeId{9}).GetMutableEdgeMutations()[{3}] = {'A', 'T'};

  dag.RecomputeCompactGenomes();
  return input_storage;
}

#pragma once

#include "larch/madag/mutation_annotated_dag.hpp"

using NodeSeqMap = std::unordered_map<NodeId, std::string>;

[[maybe_unused]] static auto MakeSampleDAG() {
  MADAGStorage input_storage{{}};
  auto dag = input_storage.View();

  dag.SetReferenceSequence("GAA");
  dag.InitializeNodes(11);

  size_t edge_id = 0;
  dag.AddEdge({edge_id++}, {0}, {10}, {0});
  dag.AddEdge({edge_id++}, {7}, {1}, {0}).GetMutableEdgeMutations()[{1}] = {'T', 'A'};
  dag.AddEdge({edge_id++}, {7}, {2}, {1}).GetMutableEdgeMutations()[{1}] = {'T', 'G'};
  dag.AddEdge({edge_id++}, {8}, {3}, {0}).GetMutableEdgeMutations()[{1}] = {'C', 'A'};
  dag.AddEdge({edge_id++}, {8}, {4}, {1}).GetMutableEdgeMutations()[{1}] = {'C', 'A'};
  dag.AddEdge({edge_id++}, {9}, {5}, {0}).GetMutableEdgeMutations()[{1}] = {'A', 'C'};
  dag.AddEdge({edge_id++}, {9}, {6}, {1}).GetMutableEdgeMutations()[{1}] = {'A', 'T'};
  dag.AddEdge({edge_id++}, {8}, {7}, {2}).GetMutableEdgeMutations()[{1}] = {'C', 'T'};
  dag.AddEdge({edge_id++}, {10}, {8}, {0}).GetMutableEdgeMutations()[{1}] = {'G', 'C'};
  dag.AddEdge({edge_id++}, {10}, {9}, {1}).GetMutableEdgeMutations()[{1}] = {'G', 'A'};

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

  dag.RecomputeCompactGenomes(true);
  dag.SampleIdsFromCG();

  return input_storage;
}

[[maybe_unused]] static auto MakeSampleDAGTopology() {
  MADAGStorage dag_storage{{}};
  auto dag = dag_storage.View();

  dag.SetReferenceSequence("GAA");

  dag.InitializeNodes(11);

  size_t edge_id = 0;
  dag.AddEdge({edge_id++}, {0}, {10}, {0});
  dag.AddEdge({edge_id++}, {7}, {1}, {0});
  dag.AddEdge({edge_id++}, {7}, {2}, {1});
  dag.AddEdge({edge_id++}, {8}, {3}, {0});
  dag.AddEdge({edge_id++}, {8}, {4}, {1});
  dag.AddEdge({edge_id++}, {9}, {5}, {0});
  dag.AddEdge({edge_id++}, {9}, {6}, {1});
  dag.AddEdge({edge_id++}, {8}, {7}, {2});
  dag.AddEdge({edge_id++}, {10}, {8}, {0});
  dag.AddEdge({edge_id++}, {10}, {9}, {1});

  dag.BuildConnections();

  return dag_storage;
}

[[maybe_unused]] static auto MakeSampleUnambiguousCompleteSequenceMap() {
  NodeSeqMap node_seq_map;
  // leaf nodes
  node_seq_map[{1}] = {"ACC"};
  node_seq_map[{2}] = {"TAG"};
  node_seq_map[{3}] = {"GGG"};
  node_seq_map[{4}] = {"ACG"};
  node_seq_map[{5}] = {"CTT"};
  node_seq_map[{6}] = {"TCC"};
  // internal nodes
  node_seq_map[{7}] = {"TGG"};
  node_seq_map[{8}] = {"GTC"};
  node_seq_map[{9}] = {"AGT"};
  node_seq_map[{10}] = {"GAA"};
  return node_seq_map;
}

[[maybe_unused]] static auto MakeSampleUnambiguousCompleteSequenceMap2() {
  NodeSeqMap node_seq_map;
  // leaf nodes
  node_seq_map[{1}] = {"ACC"};
  node_seq_map[{2}] = {"TAG"};
  node_seq_map[{3}] = {"GGG"};
  node_seq_map[{4}] = {"ACG"};
  node_seq_map[{5}] = {"CTT"};
  node_seq_map[{6}] = {"TCC"};
  // internal nodes
  node_seq_map[{7}] = {"TGG"};
  node_seq_map[{8}] = {"CTC"};
  node_seq_map[{9}] = {"AGT"};
  node_seq_map[{10}] = {"GAA"};
  return node_seq_map;
}

[[maybe_unused]] static auto MakeSampleAmbiguousCompleteSequenceMap() {
  NodeSeqMap node_seq_map = MakeSampleUnambiguousCompleteSequenceMap();
  node_seq_map[{2}] = {"TNN"};
  node_seq_map[{4}] = {"ANG"};
  return node_seq_map;
}

[[maybe_unused]] static auto MakeAmbiguousSampleDAG() {
  auto amb_dag_storage = MakeSampleDAGTopology();
  auto amb_seq_map = MakeSampleAmbiguousCompleteSequenceMap();
  auto amb_dag = amb_dag_storage.View();
  amb_dag.SetCompactGenomesFromNodeSequenceMap(amb_seq_map);
  amb_dag.RecomputeEdgeMutations();
  return amb_dag_storage;
}

[[maybe_unused]] static auto MakeUnambiguousSampleDAG() {
  auto unamb_dag_storage = MakeSampleDAGTopology();
  auto unamb_seq_map = MakeSampleUnambiguousCompleteSequenceMap();
  auto unamb_dag = unamb_dag_storage.View();
  unamb_dag.SetCompactGenomesFromNodeSequenceMap(unamb_seq_map);
  unamb_dag.RecomputeEdgeMutations();
  return unamb_dag_storage;
}

[[maybe_unused]] static auto MakeUnambiguousSampleDAG2() {
  auto unamb_dag_storage = MakeSampleDAGTopology();
  auto unamb_seq_map = MakeSampleUnambiguousCompleteSequenceMap2();
  auto unamb_dag = unamb_dag_storage.View();
  unamb_dag.SetCompactGenomesFromNodeSequenceMap(unamb_seq_map);
  unamb_dag.RecomputeEdgeMutations();
  return unamb_dag_storage;

  [[maybe_unused]] static auto MakeBigSampleDAGTopology(bool missing_edges = false) {
    MADAGStorage dag_storage{{}};
    auto dag = dag_storage.View();

    dag.SetReferenceSequence("GAA");
    dag.InitializeNodes(18);

    size_t edge_id = 0;
    dag.AddEdge({edge_id++}, {0}, {17}, {0});
    dag.AddEdge({edge_id++}, {0}, {16}, {0});
    dag.AddEdge({edge_id++}, {17}, {2}, {0});
    dag.AddEdge({edge_id++}, {17}, {15}, {1});
    dag.AddEdge({edge_id++}, {16}, {1}, {0});
    dag.AddEdge({edge_id++}, {16}, {14}, {1});
    dag.AddEdge({edge_id++}, {15}, {1}, {0});
    if (!missing_edges) dag.AddEdge({edge_id++}, {15}, {13}, {1});
    dag.AddEdge({edge_id++}, {15}, {12}, {1});
    dag.AddEdge({edge_id++}, {14}, {2}, {0});
    dag.AddEdge({edge_id++}, {14}, {13}, {1});
    if (!missing_edges) dag.AddEdge({edge_id++}, {14}, {12}, {1});
    dag.AddEdge({edge_id++}, {13}, {4}, {0});
    dag.AddEdge({edge_id++}, {13}, {11}, {1});
    dag.AddEdge({edge_id++}, {12}, {3}, {0});
    dag.AddEdge({edge_id++}, {12}, {10}, {1});
    dag.AddEdge({edge_id++}, {11}, {3}, {0});
    if (!missing_edges) dag.AddEdge({edge_id++}, {11}, {9}, {1});
    dag.AddEdge({edge_id++}, {11}, {8}, {1});
    dag.AddEdge({edge_id++}, {10}, {4}, {0});
    if (!missing_edges) dag.AddEdge({edge_id++}, {10}, {8}, {1});
    dag.AddEdge({edge_id++}, {10}, {9}, {1});
    dag.AddEdge({edge_id++}, {9}, {5}, {0});
    dag.AddEdge({edge_id++}, {9}, {6}, {1});
    dag.AddEdge({edge_id++}, {9}, {7}, {2});
    dag.AddEdge({edge_id++}, {8}, {5}, {0});
    dag.AddEdge({edge_id++}, {8}, {6}, {1});
    dag.AddEdge({edge_id++}, {8}, {7}, {2});

    dag.BuildConnections();

    return dag_storage;
  }

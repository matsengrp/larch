#pragma once

#include "larch/madag/mutation_annotated_dag.hpp"

// DAG Utilities

[[maybe_unused]] static void dag_info(MADAGStorage<>& dag_storage,
                                      bool summary_only = true) {
  auto dag = dag_storage.View();
  std::cout << "=== DAG_INFO [begin] ===" << std::endl;
  std::cout << "=== REF_SEQ: " << dag.GetReferenceSequence().size()
            << " ===" << std::endl;
  std::cout << "=== NODES: " << dag.GetNodesCount() << " ===" << std::endl;
  std::cout << "=== EDGES: " << dag.GetEdgesCount() << " ===" << std::endl;
  std::cout << "=== LEAFS: " << dag.GetLeafsCount() << " ===" << std::endl;
  if (!summary_only) {
    std::cout << "=== REF_SEQ: " << dag.GetReferenceSequence().size()
              << " ===" << std::endl;
    std::cout << "  " << dag.GetReferenceSequence() << std::endl;
    std::cout << "=== NODES: " << dag.GetNodesCount() << " ===" << std::endl;
    for (auto node : dag.GetNodes()) {
      std::cout << "node: " << node.GetId() << std::endl;
      std::cout << "  sample_id: " << node.GetSampleId().value_or("<NO_ID>")
                << std::endl;
      std::cout << "  leaf_set: " << node.GetLeafsBelow() << std::endl;
      std::cout << "  children: " << node.GetChildren() << std::endl;
      std::cout << "  clades: " << node.GetClades() << std::endl;
    }
    std::cout << "=== EDGES: " << dag.GetEdgesCount() << " ===" << std::endl;
    for (auto edge : dag.GetEdges()) {
      std::cout << "edge: " << edge.GetId() << " [" << edge.GetParent() << " -> "
                << edge.GetChild() << " | " << edge.GetClade() << "]" << std::endl;
      std::cout << "  muts: [";
      for (auto [pos, nucs] : edge.GetEdgeMutations()) {
        std::cout << "" << nucs.first.ToChar() << pos << nucs.second.ToChar() << ",";
      }
      std::cout << "]" << std::endl;
    }
  }
  std::cout << "=== DAG_INFO [end] ===" << std::endl << std::endl;
}

template <typename DAG1, typename DAG2>
bool compare_treedags(DAG1 dag1, DAG2 dag2, bool do_print = true) {
  if (dag1.GetReferenceSequence() != dag2.GetReferenceSequence()) {
    if (do_print) {
      std::cout << "DAG_MISMATCH: ref_seq" << std::endl;
    }
    return false;
  }
  if (dag1.GetNodesCount() != dag2.GetNodesCount()) {
    if (do_print) {
      std::cout << "DAG_MISMATCH: node_count" << std::endl;
    }
    return false;
  }
  if (dag1.GetEdgesCount() != dag2.GetEdgesCount()) {
    if (do_print) {
      std::cout << "DAG_MISMATCH: edge_count" << std::endl;
    }
    return false;
  }

  std::unordered_set<CompactGenome> dag1_cgs, dag2_cgs;
  for (auto node : dag1.GetNodes()) {
    dag1_cgs.emplace(node.GetCompactGenome().Copy(&node));
  }
  for (auto node : dag2.GetNodes()) {
    dag2_cgs.emplace(node.GetCompactGenome().Copy(&node));
  }
  if (dag1_cgs != dag2_cgs) {
    if (do_print) {
      std::cout << "DAG_MISMATCH: cgs" << std::endl;
    }
    return false;
  }

  std::vector<EdgeMutations> dag1_ems;
  std::vector<EdgeMutations> dag2_ems;

  for (auto edge : dag1.GetEdges()) {
    dag1_ems.emplace_back(edge.GetEdgeMutations().Copy(&edge));
  }
  for (auto edge : dag2.GetEdges()) {
    dag2_ems.emplace_back(edge.GetEdgeMutations().Copy(&edge));
  }

  if (not(dag1_ems.empty() or dag2_ems.empty())) {
    for (auto& em : dag1_ems) {
      if (std::count(dag2_ems.begin(), dag2_ems.end(), em) !=
          std::count(dag1_ems.begin(), dag1_ems.end(), em)) {
        if (do_print) {
          std::cout << "DAG_MISMATCH: ems" << std::endl;
        }
        return false;
      }
    }
  } else if (not(dag1_ems.empty() and dag2_ems.empty())) {
    if (do_print) {
      std::cout << "DAG_MISMATCH: ems empty" << std::endl;
    }
    return false;
  }
  return true;
}

// Sample DAGs

[[maybe_unused]] static auto make_sample_dag() {
  MADAGStorage<> input_storage = MADAGStorage<>::EmptyDefault();
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

[[maybe_unused]] static auto make_sample_dag_topology() {
  MADAGStorage<> dag_storage = MADAGStorage<>::EmptyDefault();
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

[[maybe_unused]] static auto make_sample_sequence_map() {
  NodeSeqMap node_seq_map;
  // leaf nodes
  node_seq_map[{1}] = {"ACC"};
  node_seq_map[{2}] = {"GTT"};
  node_seq_map[{3}] = {"AGG"};
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

[[maybe_unused]] static auto make_sample_unambiguous_sequence_map() {
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

[[maybe_unused]] static auto make_sample_ambiguous_sequence_map() {
  NodeSeqMap node_seq_map = make_sample_unambiguous_sequence_map();
  // adding ambiguity to leaf nodes.
  node_seq_map[{2}] = {"TNN"};
  node_seq_map[{4}] = {"ANG"};
  return node_seq_map;
}

[[maybe_unused]] static auto make_sample_unambiguous_sequence_map_2() {
  NodeSeqMap node_seq_map = make_sample_unambiguous_sequence_map();
  // alternate leaf nodes
  node_seq_map[{1}] = {"TGA"};  // no matches
  node_seq_map[{2}] = {"AAG"};  // pos 1 differing
  node_seq_map[{3}] = {"GCG"};  // pos 2 differing
  node_seq_map[{4}] = {"CGA"};  // shuffled
  node_seq_map[{5}] = {"CGG"};  // 2 differing
  node_seq_map[{6}] = {"TCC"};  // full match
  return node_seq_map;
}

[[maybe_unused]] static auto make_base_sample_dag() {
  auto dag_storage = make_sample_dag_topology();
  auto seq_map = make_sample_sequence_map();
  auto dag = dag_storage.View();
  dag.SetCompactGenomesFromNodeSequenceMap(seq_map);
  dag.RecomputeEdgeMutations();
  dag.RecomputeCompactGenomes(false);
  dag.SampleIdsFromCG();
  return dag_storage;
}

[[maybe_unused]] static auto make_ambiguous_sample_dag() {
  auto dag_storage = make_sample_dag_topology();
  auto seq_map = make_sample_ambiguous_sequence_map();
  auto dag = dag_storage.View();
  dag.SetCompactGenomesFromNodeSequenceMap(seq_map);
  dag.RecomputeEdgeMutations();
  dag.RecomputeCompactGenomes(false);
  dag.SampleIdsFromCG();
  return dag_storage;
}

[[maybe_unused]] static auto make_unambiguous_sample_dag() {
  auto dag_storage = make_sample_dag_topology();
  auto seq_map = make_sample_unambiguous_sequence_map();
  auto dag = dag_storage.View();
  dag.SetCompactGenomesFromNodeSequenceMap(seq_map);
  dag.RecomputeEdgeMutations();
  dag.RecomputeCompactGenomes(false);
  dag.SampleIdsFromCG();
  return dag_storage;
}

[[maybe_unused]] static auto make_unambiguous_sample_dag_2() {
  auto dag_storage = make_sample_dag_topology();
  auto seq_map = make_sample_unambiguous_sequence_map_2();
  auto dag = dag_storage.View();
  dag.SetCompactGenomesFromNodeSequenceMap(seq_map);
  dag.RecomputeEdgeMutations();
  dag.RecomputeCompactGenomes(false);
  dag.SampleIdsFromCG();
  return dag_storage;
}

[[maybe_unused]] static auto make_nonintersecting_sample_dag_topology() {
  MADAGStorage<> dag_storage = MADAGStorage<>::EmptyDefault();
  auto dag = dag_storage.View();
  dag.SetReferenceSequence("GAA");
  dag.InitializeNodes(11);

  size_t edge_id = 0;
  dag.AddEdge({edge_id++}, {0}, {10}, {0});
  dag.AddEdge({edge_id++}, {8}, {1}, {0});
  dag.AddEdge({edge_id++}, {8}, {2}, {1});
  dag.AddEdge({edge_id++}, {10}, {3}, {0});
  dag.AddEdge({edge_id++}, {7}, {4}, {0});
  dag.AddEdge({edge_id++}, {7}, {5}, {1});
  dag.AddEdge({edge_id++}, {9}, {6}, {0});
  dag.AddEdge({edge_id++}, {9}, {7}, {1});
  dag.AddEdge({edge_id++}, {10}, {8}, {1});
  dag.AddEdge({edge_id++}, {10}, {9}, {2});

  dag.BuildConnections();
  return dag_storage;
}

[[maybe_unused]] static auto make_nonintersectioning_sequence_map() {
  NodeSeqMap node_seq_map;
  // leaf nodes
  node_seq_map[{1}] = {"ACC"};
  node_seq_map[{2}] = {"TAG"};
  node_seq_map[{3}] = {"GAA"};
  node_seq_map[{4}] = {"ACG"};
  node_seq_map[{5}] = {"CTT"};
  node_seq_map[{6}] = {"TCC"};
  // internal nodes
  node_seq_map[{7}] = {"TGG"};
  node_seq_map[{8}] = {"GAA"};
  node_seq_map[{9}] = {"GAA"};
  node_seq_map[{10}] = {"GAA"};
  return node_seq_map;
}

// `missing_edges` option determines whether to build a complete or incomplete DAG.
// Complete DAG by default.
[[maybe_unused]] static auto make_big_sample_dag_topology(
    bool is_missing_edges = false) {
  MADAGStorage<> dag_storage = MADAGStorage<>::EmptyDefault();
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
  if (!is_missing_edges) dag.AddEdge({edge_id++}, {15}, {13}, {1});
  dag.AddEdge({edge_id++}, {15}, {12}, {1});
  dag.AddEdge({edge_id++}, {14}, {2}, {0});
  dag.AddEdge({edge_id++}, {14}, {13}, {1});
  if (!is_missing_edges) dag.AddEdge({edge_id++}, {14}, {12}, {1});
  dag.AddEdge({edge_id++}, {13}, {4}, {0});
  dag.AddEdge({edge_id++}, {13}, {11}, {1});
  dag.AddEdge({edge_id++}, {12}, {3}, {0});
  dag.AddEdge({edge_id++}, {12}, {10}, {1});
  dag.AddEdge({edge_id++}, {11}, {3}, {0});
  if (!is_missing_edges) dag.AddEdge({edge_id++}, {11}, {9}, {1});
  dag.AddEdge({edge_id++}, {11}, {8}, {1});
  dag.AddEdge({edge_id++}, {10}, {4}, {0});
  if (!is_missing_edges) dag.AddEdge({edge_id++}, {10}, {8}, {1});
  dag.AddEdge({edge_id++}, {10}, {9}, {1});
  dag.AddEdge({edge_id++}, {9}, {5}, {0});
  dag.AddEdge({edge_id++}, {9}, {6}, {1});
  dag.AddEdge({edge_id++}, {9}, {7}, {2});
  dag.AddEdge({edge_id++}, {8}, {5}, {0});
  dag.AddEdge({edge_id++}, {8}, {6}, {1});
  dag.AddEdge({edge_id++}, {8}, {7}, {2});
  dag.BuildConnections();

  dag.RecomputeCompactGenomes(true);
  dag.SampleIdsFromCG();
  return dag_storage;
}

// Sample DAG, but contains one addition

[[maybe_unused]] static auto make_sample_dag_with_one_unique_node_topology() {
  MADAGStorage<> dag_storage = MADAGStorage<>::EmptyDefault();
  auto dag = dag_storage.View();
  dag.SetReferenceSequence("GAA");
  dag.InitializeNodes(12);

  size_t edge_id = 0;
  dag.AddEdge({edge_id++}, {0}, {11}, {0});
  dag.AddEdge({edge_id++}, {7}, {1}, {0});
  dag.AddEdge({edge_id++}, {7}, {2}, {1});
  dag.AddEdge({edge_id++}, {8}, {3}, {0});
  dag.AddEdge({edge_id++}, {8}, {4}, {1});
  dag.AddEdge({edge_id++}, {9}, {5}, {0});
  dag.AddEdge({edge_id++}, {9}, {6}, {1});
  dag.AddEdge({edge_id++}, {10}, {8}, {0});   // New node
  dag.AddEdge({edge_id++}, {10}, {7}, {1});   // New node
  dag.AddEdge({edge_id++}, {11}, {10}, {0});  // 10 updated to 11
  dag.AddEdge({edge_id++}, {11}, {9}, {1});   // 10 updated to 11

  dag.BuildConnections();
  return dag_storage;
}

[[maybe_unused]] static auto make_sample_dag_with_one_unique_node_sequence_map() {
  NodeSeqMap node_seq_map;
  // leaf nodes
  node_seq_map[{1}] = {"ACC"};
  node_seq_map[{2}] = {"GTT"};
  node_seq_map[{3}] = {"AGG"};
  node_seq_map[{4}] = {"ACG"};
  node_seq_map[{5}] = {"CTT"};
  node_seq_map[{6}] = {"TCC"};
  // internal nodes
  node_seq_map[{7}] = {"TGG"};
  node_seq_map[{8}] = {"CTC"};
  node_seq_map[{9}] = {"AGT"};
  node_seq_map[{10}] = {"GGA"};  // New node.
  node_seq_map[{11}] = {"GAA"};  // 10 updated to 11.
  return node_seq_map;
}

[[maybe_unused]] static auto make_sample_dag_with_one_unique_node() {
  auto dag_storage = make_sample_dag_with_one_unique_node_topology();
  auto seq_map = make_sample_dag_with_one_unique_node_sequence_map();
  auto dag = dag_storage.View();
  dag.SetCompactGenomesFromNodeSequenceMap(seq_map);
  dag.RecomputeEdgeMutations();
  dag.RecomputeCompactGenomes(false);
  dag.SampleIdsFromCG();
  return dag_storage;
}

// DAG which does not intersect with Sample DAG.

[[maybe_unused]] static auto make_nonintersecting_sample_dag() {
  MADAGStorage<> input_storage = MADAGStorage<>::EmptyDefault();
  auto dag = input_storage.View();
  dag.SetReferenceSequence("GAA");
  dag.InitializeNodes(11);
  dag.AddEdge({0}, {0}, {10}, {0});
  dag.AddEdge({1}, {8}, {1}, {0}).GetMutableEdgeMutations()[{1}] = {'A', 'T'};
  dag.AddEdge({2}, {8}, {2}, {1}).GetMutableEdgeMutations()[{2}] = {'A', 'C'};
  dag.AddEdge({3}, {10}, {3}, {0}).GetMutableEdgeMutations()[{2}] = {'A', 'T'};
  dag.AddEdge({4}, {7}, {4}, {0}).GetMutableEdgeMutations()[{2}] = {'A', 'G'};
  dag.AddEdge({5}, {7}, {5}, {1}).GetMutableEdgeMutations()[{2}] = {'A', 'C'};
  dag.AddEdge({6}, {9}, {6}, {0}).GetMutableEdgeMutations()[{1}] = {'A', 'C'};
  dag.AddEdge({7}, {9}, {7}, {1});
  dag.AddEdge({8}, {10}, {8}, {1}).GetMutableEdgeMutations()[{1}] = {'G', 'A'};
  dag.AddEdge({9}, {10}, {9}, {2}).GetMutableEdgeMutations()[{1}] = {'G', 'A'};
  dag.BuildConnections();
  dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{2}] = {'A', 'C'};
  dag.Get(EdgeId{6}).GetMutableEdgeMutations()[{2}] = {'A', 'T'};
  dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{3}] = {'A', 'C'};
  dag.Get(EdgeId{2}).GetMutableEdgeMutations()[{3}] = {'A', 'G'};
  dag.Get(EdgeId{3}).GetMutableEdgeMutations()[{3}] = {'A', 'T'};
  dag.Get(EdgeId{4}).GetMutableEdgeMutations()[{3}] = {'A', 'G'};
  dag.Get(EdgeId{5}).GetMutableEdgeMutations()[{3}] = {'A', 'C'};
  dag.Get(EdgeId{6}).GetMutableEdgeMutations()[{3}] = {'A', 'T'};
  dag.RecomputeCompactGenomes(true);
  dag.SampleIdsFromCG();
  return input_storage;
}

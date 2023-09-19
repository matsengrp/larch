#include "larch/rf_distance.hpp"

#include "test_common.hpp"
#include "sample_dag.hpp"
#include "larch/subtree/subtree_weight.hpp"

static void test_zero_rf_distance() {
  auto storage = MakeSampleDAG();
  auto view = storage.View();
  Merge merge{view.GetReferenceSequence()};
  merge.AddDAG(view);
  auto dag = merge.GetResult();
  SubtreeWeight<SumRFDistance, std::decay_t<decltype(dag)>> count{dag};
  auto score = count.ComputeWeightBelow(dag.GetRoot(), RFDistance{merge});
  std::cout << "RF distance: " << score << "\n";
}

static void test_rf_on_two_identical_topologies() {
  auto dag1_storage = MakeSampleDAG();
  auto dag2_storage = MakeSampleDAG();
  auto dag1 = dag1_storage.View();
  auto dag2 = dag2_storage.View();
  dag1.RecomputeCompactGenomes();
  dag2.RecomputeCompactGenomes();

  // change the compact genomes for internal nodes of dag2.
  dag2.Get(NodeId{7}) = CompactGenome{"AAA", "GAA"};
  dag2.Get(NodeId{8}) = CompactGenome{"AAA", "GAA"};
  dag2.Get(NodeId{9}) = CompactGenome{"AAA", "GAA"};
  dag2.Get(NodeId{10}) = CompactGenome{"AAA", "GAA"};
  dag2.RecomputeEdgeMutations();

  Merge merge(dag1.GetReferenceSequence());
  merge.AddDAGs(std::vector{dag1, dag2});

  auto dag = merge.GetResult();
  SubtreeWeight<SumRFDistance, std::decay_t<decltype(dag)>> count{dag};
  auto score = count.ComputeWeightBelow(dag.GetRoot(), RFDistance{merge});
  Assert(score == 0);
}

static MADAGStorage MakeNonintersectingSampleDAG() {
  MADAGStorage input_storage{{}};
  auto dag = input_storage.View();
  dag.SetReferenceSequence("GAA");
  dag.InitializeNodes(11);
  dag.AddEdge({0}, {0}, {10}, {0});
  dag.AddEdge({1}, {8}, {1}, {0}).GetMutableEdgeMutations()[{1}] = {'A', 'T'};
  dag.AddEdge({2}, {8}, {2}, {1}).GetMutableEdgeMutations()[{2}] = {'A', 'C'};
  dag.AddEdge({3}, {10}, {3}, {0}).GetMutableEdgeMutations()[{1}] = {'G', 'T'};
  dag.AddEdge({4}, {7}, {4}, {0}).GetMutableEdgeMutations()[{2}] = {'A', 'G'};
  dag.AddEdge({5}, {7}, {5}, {1}).GetMutableEdgeMutations()[{2}] = {'A', 'C'};
  dag.AddEdge({6}, {9}, {6}, {0}).GetMutableEdgeMutations()[{1}] = {'A', 'C'};
  dag.AddEdge({7}, {9}, {7}, {1});
  dag.AddEdge({9}, {10}, {8}, {1}).GetMutableEdgeMutations()[{1}] = {'G', 'A'};
  dag.AddEdge({8}, {10}, {9}, {2}).GetMutableEdgeMutations()[{1}] = {'G', 'A'};
  dag.BuildConnections();
  dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{2}] = {'A', 'C'};
  dag.Get(EdgeId{3}).GetMutableEdgeMutations()[{2}] = {'A', 'T'};
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

static void test_rf_distance_hand_computed_example() {
  auto dag1_storage = MakeSampleDAG();
  auto dag2_storage = MakeNonintersectingSampleDAG();
  auto dag1 = dag1_storage.View();
  auto dag2 = dag2_storage.View();
  Merge merge(dag1.GetReferenceSequence());
  merge.AddDAGs(std::vector{dag1, dag2});
  auto dag = merge.GetResult();
  SubtreeWeight<SumRFDistance, std::decay_t<decltype(dag)>> count{dag};
  auto score = count.ComputeWeightBelow(dag.GetRoot(), RFDistance{merge});
  Assert(score == 2 * (dag1.GetNodesCount() + dag2.GetNodesCount()));
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_zero_rf_distance(); }, "RF distance: zero"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_rf_on_two_identical_topologies(); },
              "RF distance: identical topologies"});

[[maybe_unused]] static const auto test_added2 = add_test(
    {[] { test_rf_distance_hand_computed_example(); }, "RF distance: hand computed"});

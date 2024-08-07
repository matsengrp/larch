#include "test_common.hpp"
#include "test_common_dag.hpp"
#include "larch/dag_loader.hpp"

using Node = MutableMADAG::NodeView;
using Edge = MutableMADAG::EdgeView;

[[maybe_unused]] std::set<std::tuple<NodeId, NodeId, CladeIdx>> dag_make_node_pair_map(
    MADAGStorage<>& dag_storage) {
  auto dag = dag_storage.View();
  std::set<std::tuple<NodeId, NodeId, CladeIdx>> node_pairs;
  for (auto edge : dag.GetEdges()) {
    node_pairs.insert({edge.GetParent(), edge.GetChild(), edge.GetClade()});
  }
  return node_pairs;
}

[[maybe_unused]] bool dag_compare_topologies(MADAGStorage<>& lhs_storage,
                                             MADAGStorage<>& rhs_storage) {
  auto lhs_node_map = dag_make_node_pair_map(lhs_storage);
  auto rhs_node_map = dag_make_node_pair_map(rhs_storage);
  bool compare_node_maps = (lhs_node_map == rhs_node_map);
  return compare_node_maps;
}

[[maybe_unused]] void dag_make_complete_rootsplits(MADAGStorage<>& dag_storage) {
  auto dag = dag_storage.View();
  auto clade_union_map = dag.BuildCladeUnionMap();
  size_t leaf_count = 0;
  std::set<NodeId>* rootsplits = nullptr;
  for (auto& [clade_union, node_ids] : clade_union_map) {
    if (clade_union.size() > leaf_count) {
      rootsplits = &node_ids;
    }
    leaf_count = std::max(leaf_count, clade_union.size());
  }
  for (auto node_id : *rootsplits) {
    dag.AppendEdge(dag.GetRoot().GetId(), node_id, {0});
  }
}

[[maybe_unused]] void test_dag_completion_single(MADAGStorage<>& dag_storage,
                                                 MADAGStorage<>& dag_storage_truth) {
  auto dag = dag_storage.View();
  dag.GetRoot().CalculateLeafsBelow();
  // dag_info(dag_storage);
  dag.ClearConnections();
  // dag_info(dag_storage);
  TestAssert((not dag_compare_topologies(dag_storage, dag_storage_truth)) &&
             "DAGs are incorrectly equal after removing edges.");

  // dag_make_complete_rootsplits(dag_storage);
  dag.MakeComplete();
  // dag_info(dag_storage);
  TestAssert(dag_compare_topologies(dag_storage, dag_storage_truth) &&
             "DAGs are not equal after recompleting the DAG.");
}

[[maybe_unused]] static void test_sample_dag_completion() {
  auto dag_storage = make_sample_dag_topology();
  auto dag_storage_truth = make_sample_dag_topology();
  TestAssert(dag_compare_topologies(dag_storage, dag_storage_truth) &&
             "DAGs are not equal before altering.");
  test_dag_completion_single(dag_storage, dag_storage_truth);
}

[[maybe_unused]] static void test_big_sample_dag_completion_without_missing_edges() {
  auto big_dag_storage_complete = make_big_sample_dag_topology(false);
  auto big_dag_storage_truth = make_big_sample_dag_topology();
  TestAssert(dag_compare_topologies(big_dag_storage_complete, big_dag_storage_truth) &&
             "DAGs are not equal before altering.");
  test_dag_completion_single(big_dag_storage_complete, big_dag_storage_truth);
}

[[maybe_unused]] static void test_big_sample_dag_completion_with_missing_edges() {
  auto big_dag_storage_incomplete = make_big_sample_dag_topology(true);
  auto big_dag_storage_truth = make_big_sample_dag_topology();
  TestAssert(
      (not dag_compare_topologies(big_dag_storage_incomplete, big_dag_storage_truth)) &&
      "DAGs are incorrectly equal before altering.");
  test_dag_completion_single(big_dag_storage_incomplete, big_dag_storage_truth);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_sample_dag_completion(); }, "DAG Completion: Sample DAG"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_big_sample_dag_completion_without_missing_edges(); },
              "DAG Completion: Big Sample DAG (without missing edges)"});

[[maybe_unused]] static const auto test_added2 =
    add_test({[] { test_big_sample_dag_completion_with_missing_edges(); },
              "DAG Completion: Big Sample DAG (with missing edges)"});

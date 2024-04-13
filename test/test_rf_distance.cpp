#include "larch/rf_distance.hpp"

#include "test_common.hpp"
#include "sample_dag.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/dag_loader.hpp"

enum class RFDistanceType { Min, MinSum, Max, MaxSum };

[[maybe_unused]] static auto get_rf_distance(
    const Merge& comp_merge1, const Merge& ref_merge2,
    RFDistanceType rf_dist_type = RFDistanceType::MinSum, bool print_info = true) {
  // comp_merge1 is the DAG we compute the weights for (summing distances to ref_merge2)
  // ref_merge2 is the reference DAG.
  auto dag1 = comp_merge1.GetResult();
  if (print_info) {
    std::cout << "dag1 (compute dag) address: " << &dag1.GetStorage() << "\n";
    std::cout << "dag2 (reference dag) address: "
              << &ref_merge2.GetResult().GetStorage() << "\n";
  }
  ArbitraryInt shift_sum, result;

  // compute min rf_distance
  if (rf_dist_type == RFDistanceType::Min) {
    SubtreeWeight<RFDistance, MergeDAG> count{dag1};
    RFDistance weight_ops{ref_merge2, comp_merge1};
    shift_sum = weight_ops.GetOps().GetShiftSum();
    result =
        count.ComputeWeightBelow(dag1.GetRoot(), std::move(weight_ops)) + shift_sum;
  }
  // compute minsum rf_distance
  else if (rf_dist_type == RFDistanceType::MinSum) {
    SubtreeWeight<SumRFDistance, MergeDAG> count{dag1};
    SumRFDistance weight_ops{ref_merge2, comp_merge1};
    shift_sum = weight_ops.GetOps().GetShiftSum();
    result =
        count.ComputeWeightBelow(dag1.GetRoot(), std::move(weight_ops)) + shift_sum;
  }
  // compute max rf_distance
  else if (rf_dist_type == RFDistanceType::Max) {
    SubtreeWeight<MaxRFDistance, MergeDAG> count{dag1};
    MaxRFDistance weight_ops{ref_merge2, comp_merge1};
    shift_sum = weight_ops.GetOps().GetShiftSum();
    result =
        count.ComputeWeightBelow(dag1.GetRoot(), std::move(weight_ops)) + shift_sum;
  }
  // compute maxsum rf_distance
  else if (rf_dist_type == RFDistanceType::MaxSum) {
    SubtreeWeight<MaxSumRFDistance, MergeDAG> count{dag1};
    MaxSumRFDistance weight_ops{ref_merge2, comp_merge1};
    shift_sum = weight_ops.GetOps().GetShiftSum();
    result =
        count.ComputeWeightBelow(dag1.GetRoot(), std::move(weight_ops)) + shift_sum;
  } else {
    std::cerr << "ERROR: Invalid RFDistanceType" << std::endl;
    Assert(false);
  }

  // make sure shift sum is correct
  if (ref_merge2.GetResult().IsTree()) {
    Assert(shift_sum == ref_merge2.GetResult().GetNodesCount() - 1);
  }
  if (print_info) {
    std::cout << "shift_sum: " << shift_sum << "\n";
    std::cout << "rf_distance: " << result << "\n";
  }
  return result;
}

static void test_zero_rf_distance() {
  auto storage = make_sample_dag();
  auto view = storage.View();
  Merge merge{view.GetReferenceSequence()};
  merge.AddDAG(view);
  Assert(get_rf_distance(merge, merge) == 0);
}

static void test_rf_on_two_identical_topologies() {
  auto dag1_storage = make_sample_dag();
  auto dag2_storage = make_sample_dag();
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

  Merge merge1(dag1.GetReferenceSequence());
  merge1.AddDAGs(std::vector{dag1});
  Merge merge2(dag1.GetReferenceSequence());
  merge2.AddDAGs(std::vector{dag1});
  merge2.AddDAGs(std::vector{dag2});

  Assert(get_rf_distance(merge1, merge2) == 0);
}

static void test_rf_two_distinct_topologies_single_merge() {
  auto dag1_storage = make_sample_dag();
  auto dag2_storage = make_nonintersecting_sample_dag();
  auto dag1 = dag1_storage.View();
  auto dag2 = dag2_storage.View();

  Merge merge(dag1.GetReferenceSequence());
  merge.AddDAGs(std::vector{dag1, dag2});
  auto dist = get_rf_distance(merge, merge);
  auto truedist = (dag1.GetEdgesCount() + dag2.GetEdgesCount() -
                   dag1.GetLeafs().size() - dag2.GetLeafs().size() - 2);
  if (dist != truedist) {
    std::cout << "expected distance of " << truedist << " but computed distance was "
              << dist << "\n";
    Assert(false)
  }
}

static void test_rf_distance_hand_computed_example() {
  auto dag1_storage = make_sample_dag();
  auto dag2_storage = make_nonintersecting_sample_dag();
  auto dag1 = dag1_storage.View();
  auto dag2 = dag2_storage.View();

  Merge merge1(dag1.GetReferenceSequence());
  merge1.AddDAGs(std::vector{dag1});
  Merge merge2(dag1.GetReferenceSequence());
  merge2.AddDAGs(std::vector{dag1, dag2});
  auto true_dist = dag1.GetEdgesCount() + dag2.GetEdgesCount() -
                   dag1.GetLeafs().size() - dag2.GetLeafs().size() - 2;
  Assert(get_rf_distance(merge1, merge2) == true_dist);

  Merge merge(dag1.GetReferenceSequence());
  merge.AddDAGs(std::vector{dag1, dag2});
  Assert(get_rf_distance(merge1, merge2) == get_rf_distance(merge, merge));
}

static void test_rf_distance_different_weight_ops() {
  auto dag1_storage = make_sample_dag();
  auto dag2_storage = make_nonintersecting_sample_dag();
  auto dag3_storage = make_sample_dag_with_one_unique_node();
  auto dag1 = dag1_storage.View();
  auto dag2 = dag2_storage.View();
  auto dag3 = dag3_storage.View();

  Merge merge1(dag1.GetReferenceSequence());
  merge1.AddDAGs(std::vector{dag1});
  Merge merge2(dag2.GetReferenceSequence());
  merge2.AddDAGs(std::vector{dag2});
  Merge merge3(dag1.GetReferenceSequence());
  merge3.AddDAGs(std::vector{dag3});

  Merge merge1_2(dag1.GetReferenceSequence());
  merge1_2.AddDAGs(std::vector{dag1, dag2});
  Merge merge1_3(dag1.GetReferenceSequence());
  merge1_3.AddDAGs(std::vector{dag1, dag3});
  Merge merge2_3(dag1.GetReferenceSequence());
  merge2_3.AddDAGs(std::vector{dag2, dag3});
  Merge merge1_2_3(dag1.GetReferenceSequence());
  merge1_2_3.AddDAGs(std::vector{dag1, dag2, dag3});

  // Hand-computed RF-distances for single MATs.
  std::map<std::pair<int, int>, ArbitraryInt> true_dist_map;
  true_dist_map[{1, 1}] = true_dist_map[{2, 2}] = true_dist_map[{3, 3}] = 0;
  true_dist_map[{1, 2}] = true_dist_map[{2, 1}] = 6;
  true_dist_map[{1, 3}] = true_dist_map[{3, 1}] = 1;
  true_dist_map[{2, 3}] = true_dist_map[{3, 2}] = 7;

  // Use hand-computed distances to find rf-distance types on merged MADAGs.
  auto compute_true_dist = [&](std::vector<int> compute_ids, std::vector<int> ref_ids,
                               RFDistanceType rf_dist_type, bool do_print = false) {
    ArbitraryInt total = 0;
    std::vector<ArbitraryInt> vec;

    for (auto compute_id : compute_ids) {
      std::vector<ArbitraryInt> subvec;
      for (auto ref_id : ref_ids) {
        subvec.push_back(true_dist_map[{compute_id, ref_id}]);
      }
      if (rf_dist_type == RFDistanceType::Min) {
        vec.push_back(*std::min_element(subvec.begin(), subvec.end()));
      } else if (rf_dist_type == RFDistanceType::Max) {
        vec.push_back(*std::max_element(subvec.begin(), subvec.end()));
      } else {
        vec.push_back(std::accumulate(subvec.begin(), subvec.end(), ArbitraryInt{0}));
      }
    }

    if (rf_dist_type == RFDistanceType::Min or rf_dist_type == RFDistanceType::MinSum) {
      total = *std::min_element(vec.begin(), vec.end());
    } else {
      total = *std::max_element(vec.begin(), vec.end());
    }

    if (do_print) {
      std::cout << "true_rf_distance: " << total << std::endl;
    }
    return total;
  };

  // Compute RF-Distance of single MATs (which will have the same result for all
  // methods).
  ArbitraryInt dist;
  for (auto rf_dist_type : {RFDistanceType::Min, RFDistanceType::Max,
                            RFDistanceType::MinSum, RFDistanceType::MaxSum}) {
    bool do_print = true;
    dist = get_rf_distance(merge1, merge2, rf_dist_type, do_print);
    Assert(dist == compute_true_dist({1}, {2}, rf_dist_type, do_print));
    dist = get_rf_distance(merge2, merge1, rf_dist_type, do_print);
    Assert(dist == compute_true_dist({2}, {1}, rf_dist_type, do_print));
    dist = get_rf_distance(merge1, merge3, rf_dist_type, do_print);
    Assert(dist == compute_true_dist({1}, {3}, rf_dist_type, do_print));
    dist = get_rf_distance(merge3, merge1, rf_dist_type, do_print);
    Assert(dist == compute_true_dist({3}, {1}, rf_dist_type, do_print));
    dist = get_rf_distance(merge2, merge3, rf_dist_type, do_print);
    Assert(dist == compute_true_dist({2}, {3}, rf_dist_type, do_print));
    dist = get_rf_distance(merge3, merge2, rf_dist_type, do_print);
    Assert(dist == compute_true_dist({3}, {2}, rf_dist_type, do_print));
  }

  // Compute RF-Distance of MATs to two-MAT MADAGs (sum methods will have the same
  // result as non-sum).
  {
    bool do_print = true;
    // MADAG containing itself.
    dist = get_rf_distance(merge1_2, merge1, RFDistanceType::Min, do_print);
    Assert(dist == compute_true_dist({1, 2}, {1}, RFDistanceType::Min, do_print));
    dist = get_rf_distance(merge1_2, merge1, RFDistanceType::Max, do_print);
    Assert(dist == compute_true_dist({1, 2}, {1}, RFDistanceType::Max, do_print));
    dist = get_rf_distance(merge1_2, merge1, RFDistanceType::MinSum, do_print);
    Assert(dist == compute_true_dist({1, 2}, {1}, RFDistanceType::MinSum, do_print));
    dist = get_rf_distance(merge1_2, merge1, RFDistanceType::MaxSum, do_print);
    Assert(dist == compute_true_dist({1, 2}, {1}, RFDistanceType::MaxSum, do_print));

    // MADAG containing two differing MATs.
    dist = get_rf_distance(merge2_3, merge1, RFDistanceType::Min, do_print);
    Assert(dist == compute_true_dist({2, 3}, {1}, RFDistanceType::Min, do_print));
    dist = get_rf_distance(merge2_3, merge1, RFDistanceType::Max, do_print);
    Assert(dist == compute_true_dist({2, 3}, {1}, RFDistanceType::Max, do_print));
    dist = get_rf_distance(merge2_3, merge1, RFDistanceType::MinSum, do_print);
    Assert(dist == compute_true_dist({2, 3}, {1}, RFDistanceType::MinSum, do_print));
    dist = get_rf_distance(merge2_3, merge1, RFDistanceType::MaxSum, do_print);
    Assert(dist == compute_true_dist({2, 3}, {1}, RFDistanceType::MaxSum, do_print));
  }

  // Compute RF-Distance of MATs to three-MAT MADAGs.
  {
    bool do_print = true;
    dist = get_rf_distance(merge1_2_3, merge1, RFDistanceType::Min, do_print);
    Assert(dist == compute_true_dist({1, 2, 3}, {1}, RFDistanceType::Min, do_print));
    dist = get_rf_distance(merge1_2_3, merge1, RFDistanceType::Max, do_print);
    Assert(dist == compute_true_dist({1, 2, 3}, {1}, RFDistanceType::Max, do_print));
    dist = get_rf_distance(merge1_2_3, merge1, RFDistanceType::MinSum, do_print);
    Assert(dist == compute_true_dist({1, 2, 3}, {1}, RFDistanceType::MinSum, do_print));
    dist = get_rf_distance(merge1_2_3, merge1, RFDistanceType::MaxSum, do_print);
    Assert(dist == compute_true_dist({1, 2, 3}, {1}, RFDistanceType::MaxSum, do_print));
  }

  // Compute sum RF-Distances over MADAGs.
  {
    bool do_print = true;
    dist = get_rf_distance(merge1_2_3, merge1_2, RFDistanceType::MinSum, do_print);
    Assert(dist ==
           compute_true_dist({1, 2, 3}, {1, 2}, RFDistanceType::MinSum, do_print));
    dist = get_rf_distance(merge1_2_3, merge1_2, RFDistanceType::MaxSum, do_print);
    Assert(dist ==
           compute_true_dist({1, 2, 3}, {1, 2}, RFDistanceType::MaxSum, do_print));
    dist = get_rf_distance(merge1_2, merge1_2_3, RFDistanceType::MinSum, do_print);
    Assert(dist ==
           compute_true_dist({1, 2}, {1, 2, 3}, RFDistanceType::MinSum, do_print));
    dist = get_rf_distance(merge1_2, merge1_2_3, RFDistanceType::MaxSum, do_print);
    Assert(dist ==
           compute_true_dist({1, 2}, {1, 2, 3}, RFDistanceType::MaxSum, do_print));
    dist = get_rf_distance(merge1_2_3, merge1_2_3, RFDistanceType::MinSum, do_print);
    Assert(dist ==
           compute_true_dist({1, 2, 3}, {1, 2, 3}, RFDistanceType::MinSum, do_print));
    dist = get_rf_distance(merge1_2_3, merge1_2_3, RFDistanceType::MaxSum, do_print);
    Assert(dist ==
           compute_true_dist({1, 2, 3}, {1, 2, 3}, RFDistanceType::MaxSum, do_print));
  }
}

static void test_rf_distance_different_weight_ops() {
  auto dag0_storage = make_base_sample_dag();
  auto dag1_storage = make_sample_dag();
  auto dag2_storage = make_nonintersecting_sample_dag();
  auto dag3_storage = make_sample_dag_with_one_unique_node();
  auto dag1 = dag1_storage.View();
  auto dag2 = dag2_storage.View();
  auto dag3 = dag3_storage.View();

  Merge merge1(dag1.GetReferenceSequence());
  merge1.AddDAGs(std::vector{dag1});
  Merge merge2(dag2.GetReferenceSequence());
  merge2.AddDAGs(std::vector{dag2});
  Merge merge3(dag1.GetReferenceSequence());
  merge3.AddDAGs(std::vector{dag3});

  Merge merge1_2(dag1.GetReferenceSequence());
  merge1_2.AddDAGs(std::vector{dag1, dag2});
  Merge merge1_3(dag1.GetReferenceSequence());
  merge1_3.AddDAGs(std::vector{dag1, dag3});
  Merge merge2_3(dag1.GetReferenceSequence());
  merge2_3.AddDAGs(std::vector{dag2, dag3});
  Merge merge1_2_3(dag1.GetReferenceSequence());
  merge1_2_3.AddDAGs(std::vector{dag1, dag2, dag3});

  // Hand-computed RF-distances for single MATs.
  std::map<std::pair<int, int>, ArbitraryInt> true_dist_map;
  true_dist_map[{1, 1}] = true_dist_map[{2, 2}] = true_dist_map[{3, 3}] = 0;
  true_dist_map[{1, 2}] = true_dist_map[{2, 1}] = 6;
  true_dist_map[{1, 3}] = true_dist_map[{3, 1}] = 1;
  true_dist_map[{2, 3}] = true_dist_map[{3, 2}] = 7;

  // Use hand-computed distances to find rf-distance types on merged MADAGs.
  auto compute_true_dist = [&](std::vector<int> compute_ids, std::vector<int> ref_ids,
                               RFDistanceType rf_dist_type, bool do_print = false) {
    ArbitraryInt total = 0;
    std::vector<ArbitraryInt> vec;

    for (auto compute_id : compute_ids) {
      std::vector<ArbitraryInt> subvec;
      for (auto ref_id : ref_ids) {
        subvec.push_back(true_dist_map[{compute_id, ref_id}]);
      }
      if (rf_dist_type == RFDistanceType::Min) {
        vec.push_back(*std::min_element(subvec.begin(), subvec.end()));
      } else if (rf_dist_type == RFDistanceType::Max) {
        vec.push_back(*std::max_element(subvec.begin(), subvec.end()));
      } else {
        vec.push_back(std::accumulate(subvec.begin(), subvec.end(), ArbitraryInt{0}));
      }
    }

    if (rf_dist_type == RFDistanceType::Min or rf_dist_type == RFDistanceType::MinSum) {
      total = *std::min_element(vec.begin(), vec.end());
    } else {
      total = *std::max_element(vec.begin(), vec.end());
    }

    if (do_print) {
      std::cout << "true_rf_distance: " << total << std::endl;
    }
    return total;
  };

  // Compute RF-Distance of single MATs (which will have the same result for all
  // methods).
  ArbitraryInt dist;
  for (auto rf_dist_type : {RFDistanceType::Min, RFDistanceType::Max,
                            RFDistanceType::MinSum, RFDistanceType::MaxSum}) {
    bool do_print = true;
    dist = get_rf_distance(merge1, merge2, rf_dist_type, do_print);
    Assert(dist == compute_true_dist({1}, {2}, rf_dist_type, do_print));
    dist = get_rf_distance(merge2, merge1, rf_dist_type, do_print);
    Assert(dist == compute_true_dist({2}, {1}, rf_dist_type, do_print));
    dist = get_rf_distance(merge1, merge3, rf_dist_type, do_print);
    Assert(dist == compute_true_dist({1}, {3}, rf_dist_type, do_print));
    dist = get_rf_distance(merge3, merge1, rf_dist_type, do_print);
    Assert(dist == compute_true_dist({3}, {1}, rf_dist_type, do_print));
    dist = get_rf_distance(merge2, merge3, rf_dist_type, do_print);
    Assert(dist == compute_true_dist({2}, {3}, rf_dist_type, do_print));
    dist = get_rf_distance(merge3, merge2, rf_dist_type, do_print);
    Assert(dist == compute_true_dist({3}, {2}, rf_dist_type, do_print));
  }

  // Compute RF-Distance of MATs to two-MAT MADAGs (sum methods will have the same
  // result as non-sum).
  {
    bool do_print = true;
    // MADAG containing itself.
    dist = get_rf_distance(merge1_2, merge1, RFDistanceType::Min, do_print);
    Assert(dist == compute_true_dist({1, 2}, {1}, RFDistanceType::Min, do_print));
    dist = get_rf_distance(merge1_2, merge1, RFDistanceType::Max, do_print);
    Assert(dist == compute_true_dist({1, 2}, {1}, RFDistanceType::Max, do_print));
    dist = get_rf_distance(merge1_2, merge1, RFDistanceType::MinSum, do_print);
    Assert(dist == compute_true_dist({1, 2}, {1}, RFDistanceType::MinSum, do_print));
    dist = get_rf_distance(merge1_2, merge1, RFDistanceType::MaxSum, do_print);
    Assert(dist == compute_true_dist({1, 2}, {1}, RFDistanceType::MaxSum, do_print));

    // MADAG containing two differing MATs.
    dist = get_rf_distance(merge2_3, merge1, RFDistanceType::Min, do_print);
    Assert(dist == compute_true_dist({2, 3}, {1}, RFDistanceType::Min, do_print));
    dist = get_rf_distance(merge2_3, merge1, RFDistanceType::Max, do_print);
    Assert(dist == compute_true_dist({2, 3}, {1}, RFDistanceType::Max, do_print));
    dist = get_rf_distance(merge2_3, merge1, RFDistanceType::MinSum, do_print);
    Assert(dist == compute_true_dist({2, 3}, {1}, RFDistanceType::MinSum, do_print));
    dist = get_rf_distance(merge2_3, merge1, RFDistanceType::MaxSum, do_print);
    Assert(dist == compute_true_dist({2, 3}, {1}, RFDistanceType::MaxSum, do_print));
  }

  // Compute RF-Distance of MATs to three-MAT MADAGs.
  {
    bool do_print = true;
    dist = get_rf_distance(merge1_2_3, merge1, RFDistanceType::Min, do_print);
    Assert(dist == compute_true_dist({1, 2, 3}, {1}, RFDistanceType::Min, do_print));
    dist = get_rf_distance(merge1_2_3, merge1, RFDistanceType::Max, do_print);
    Assert(dist == compute_true_dist({1, 2, 3}, {1}, RFDistanceType::Max, do_print));
    dist = get_rf_distance(merge1_2_3, merge1, RFDistanceType::MinSum, do_print);
    Assert(dist == compute_true_dist({1, 2, 3}, {1}, RFDistanceType::MinSum, do_print));
    dist = get_rf_distance(merge1_2_3, merge1, RFDistanceType::MaxSum, do_print);
    Assert(dist == compute_true_dist({1, 2, 3}, {1}, RFDistanceType::MaxSum, do_print));
  }

  // Compute sum RF-Distances over MADAGs.
  {
    bool do_print = true;
    dist = get_rf_distance(merge1_2_3, merge1_2, RFDistanceType::MinSum, do_print);
    Assert(dist ==
           compute_true_dist({1, 2, 3}, {1, 2}, RFDistanceType::MinSum, do_print));
    dist = get_rf_distance(merge1_2_3, merge1_2, RFDistanceType::MaxSum, do_print);
    Assert(dist ==
           compute_true_dist({1, 2, 3}, {1, 2}, RFDistanceType::MaxSum, do_print));
    dist = get_rf_distance(merge1_2, merge1_2_3, RFDistanceType::MinSum, do_print);
    Assert(dist ==
           compute_true_dist({1, 2}, {1, 2, 3}, RFDistanceType::MinSum, do_print));
    dist = get_rf_distance(merge1_2, merge1_2_3, RFDistanceType::MaxSum, do_print);
    Assert(dist ==
           compute_true_dist({1, 2}, {1, 2, 3}, RFDistanceType::MaxSum, do_print));
    dist = get_rf_distance(merge1_2_3, merge1_2_3, RFDistanceType::MinSum, do_print);
    Assert(dist ==
           compute_true_dist({1, 2, 3}, {1, 2, 3}, RFDistanceType::MinSum, do_print));
    dist = get_rf_distance(merge1_2_3, merge1_2_3, RFDistanceType::MaxSum, do_print);
    Assert(dist ==
           compute_true_dist({1, 2, 3}, {1, 2, 3}, RFDistanceType::MaxSum, do_print));
  }
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_zero_rf_distance(); }, "RF distance: zero"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_rf_on_two_identical_topologies(); },
              "RF distance: identical topologies"});

[[maybe_unused]] static const auto test_added2 = add_test(
    {[] { test_rf_distance_hand_computed_example(); }, "RF distance: hand computed"});

[[maybe_unused]] static const auto test_added3 =
    add_test({[] { test_rf_two_distinct_topologies_single_merge(); },
              "RF distance: distinct topologies, one merge"});

[[maybe_unused]] static const auto test_added4 =
    add_test({[] { test_rf_distance_different_weight_ops(); },
              "RF distance: using different weight ops (Min, MinSum, Max, MaxSum)"});

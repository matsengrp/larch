#include "larch/rf_distance.hpp"

#include "test_common.hpp"
#include "sample_dag.hpp"
#include "larch/subtree/subtree_weight.hpp"

static void test_rf_distance() {
  auto storage = MakeSampleDAG();
  auto view = storage.View();
  Merge merge{view.GetReferenceSequence()};
  merge.AddDAG(view);
  auto dag = merge.GetResult();
  SubtreeWeight<SumRFDistance, std::decay_t<decltype(dag)>> count{dag};
  auto score = count.ComputeWeightBelow(dag.GetRoot(), RFDistance{merge});
  std::cout << "RF distance: " << score << "\n";
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_rf_distance(); }, "RF distance"});

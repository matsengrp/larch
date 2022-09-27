#include <string_view>

#include "test_common.hpp"

#include "dag_loader.hpp"
#include "merge.hpp"

static void test_clade_idx(std::string_view path) {
  MADAG dag = LoadDAGFromProtobuf(path);
  Merge merge{dag.GetReferenceSequence()};
  merge.AddDAGs({dag});
  MADAG& merged = merge.GetResult();
  merged.GetEdgeMutations() = merge.ComputeResultEdgeMutations();
  merged.GetCompactGenomes() =
      merged.ComputeCompactGenomes(merged.GetReferenceSequence());
}

[[maybe_unused]] static const auto test_added_0 = add_test(
    {[] { test_clade_idx("data/testcase/full_dag.pb.gz"); }, "CladeIdx: testcase"});

[[maybe_unused]] static const auto test_added_1 = add_test(
    {[] { test_clade_idx("data/testcase1/full_dag.pb.gz"); }, "CladeIdx: testcase1"});

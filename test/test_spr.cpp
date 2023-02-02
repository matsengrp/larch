#include "larch/madag/compact_genome.hpp"

#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/madag/mutation_annotated_dag.hpp"

#include "larch/spr/lca.hpp"
#include "larch/spr/spr_view.hpp"

[[maybe_unused]] static void test_lca(std::string_view path) {
  MADAGStorage dag_storage = LoadDAGFromProtobuf(path);

  auto [lca, path0, path1] =
      FindLCA(dag_storage.View().Get(NodeId{30}), dag_storage.View().Get(NodeId{47}));
  Assert(lca.value == 51);
}

[[maybe_unused]] static void test_overlay(std::string_view path) {
  MADAGStorage dag_storage = LoadDAGFromProtobuf(path);
  OverlayDAGStorage overlaid{std::move(dag_storage)};

  auto node = overlaid.View().Get(NodeId{15});
  Assert(not node.IsOverlaid());
  node.Overlay();
  Assert(node.IsOverlaid());
}

[[maybe_unused]] static const auto test_added0 = add_test(
    {[] { test_lca("data/test_5_trees/tree_0.pb.gz"); }, "SPR: LCA test_5_trees"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_overlay("data/test_5_trees/tree_0.pb.gz"); },
              "SPR: overlay test_5_trees"});

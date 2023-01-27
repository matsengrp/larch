#include "larch/madag/compact_genome.hpp"

#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/madag/mutation_annotated_dag.hpp"

#include "larch/spr/lca.hpp"
#include "larch/spr/spr_view.hpp"

[[maybe_unused]] static void test_spr(std::string_view path) {
  MADAGStorage dag_storage = LoadDAGFromProtobuf(path);

  dag_storage.View().GetNodes();

  auto [lca, path0, path1] =
      FindLCA(dag_storage.View().Get(NodeId{30}), dag_storage.View().Get(NodeId{47}));
  Assert(lca.value == 51);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_spr("data/test_5_trees/tree_0.pb.gz"); }, "SPR: test_5_trees"});

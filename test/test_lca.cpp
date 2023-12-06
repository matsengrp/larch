#include "test_common.hpp"
#include "larch/spr/lca.hpp"
#include "sample_dag.hpp"

template <typename DAG>
static void AssertLCA(DAG dag, size_t n0, size_t n1, size_t lca) {
  Assert(FindLCA(dag.Get(NodeId{n0}), dag.Get(NodeId{n1})).lca.value == lca);
  Assert(FindLCA(dag.Get(NodeId{n1}), dag.Get(NodeId{n0})).lca.value == lca);
}

static void test_lca() {
  auto storage = make_sample_dag();
  auto dag = storage.View();
  AssertLCA(dag, 1, 2, 7);
  AssertLCA(dag, 1, 3, 8);
  AssertLCA(dag, 2, 4, 8);
  AssertLCA(dag, 1, 5, 10);
  AssertLCA(dag, 1, 7, 7);
  AssertLCA(dag, 1, 8, 8);
  AssertLCA(dag, 5, 10, 10);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[]() { test_lca(); }, "LCA"});

#include "test_common.hpp"
#include "larch/dag/dag.hpp"

static void test_mapped_nodes() {}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_mapped_nodes(); }, "Mapped nodes"});

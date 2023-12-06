#include <cstdlib>

#include "test_common.hpp"

static void test_larch_usher(const std::string& args) {
  std::string command =
      "./larch-usher -i data/seedtree/seedtree.pb.gz -r data/seedtree/refseq.txt.gz -o "
      "test_larch_usher_output.pb -c 2 ";
  command += args;
  TestAssert(0 == std::system(command.c_str()));
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_larch_usher(""); }, "Larch-Usher: baseline", {"slow"}});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_larch_usher("--sample-any-tree"); },
              "Larch-Usher: --sample-any-tree",
              {"slow"}});

[[maybe_unused]] static const auto test_added2 =
    add_test({[] { test_larch_usher("--switch-subtrees 2"); },
              "Larch-Usher: --switch-subtrees 2",
              {"slow"}});

[[maybe_unused]] static const auto test_added4 =
    add_test({[] { test_larch_usher("--move-coeff-nodes 0"); },
              "Larch-Usher: --move-coeff-nodes 0",
              {"slow"}});

[[maybe_unused]] static const auto test_added5 =
    add_test({[] { test_larch_usher("--callback-option best-moves"); },
              "Larch-Usher: --callback-option best-moves",
              {"slow"}});

[[maybe_unused]] static const auto test_added6 =
    add_test({[] { test_larch_usher("--callback-option best-moves-treebased"); },
              "Larch-Usher: --callback-option best-moves-treebased",
              {"slow"}});

#include <cstdlib>

#include "test_common.hpp"

static void test_larch_usher(const std::string& args) {
  std::string input_dag_path = "data/seedtree/seedtree.pb.gz";
  std::string refseq_path = "data/seedtree/refseq.txt.gz";
  std::string output_dag_path = test_output_folder + "/test_larch_usher_output.pb";
  int iter = 2;

  auto [command, result] =
      run_larch_usher(input_dag_path, output_dag_path, refseq_path, iter, args);
  TestAssert(0 == result);
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

static void test_larch_usher_merged_dag() {
  std::string output_dag_path = test_output_folder + "/opt_dag.pb";
  std::string command = "./bin/larch-usher -i data/larch_merged_dag.pb -o " +
                        output_dag_path +
                        " -c 1 -s 0 --max-subtree-clade-size 2000 --trim --quiet";
  std::cout << ">COMMAND_EXECUTE: \"" << command << "\"" << std::endl;
  int result = std::system(command.c_str());
  TestAssert(0 == result);
}

[[maybe_unused]] static const auto test_added7 = add_test(
    {test_larch_usher_merged_dag, "Larch-Usher: merged dag with trim", {"slow"}});

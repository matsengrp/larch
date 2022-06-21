#include "usher_optimize.hpp"

#include "test_common.hpp"

#include "dag_loader.hpp"

[[maybe_unused]] static void test_usher_optimize(std::string_view path) {
  InitUsherMPI(0, nullptr);

  Mutation_Annotated_Tree::Tree tree;
  Mutation_Annotated_Tree::load_mutation_annotated_tree(std::string{path}, tree);

  size_t score_before = tree.get_parsimony_score();
  UsherOptimize(tree);
  size_t score_after = tree.get_parsimony_score();

  assert_true(score_after >= score_before, "Parsimony score");
}

[[maybe_unused]] static const auto test_added = add_test(
    {[] { test_usher_optimize("data/20D_from_fasta/1final-tree-1.nh1.pb.gz"); },
     "Usher optimize"});

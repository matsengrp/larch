#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <string_view>
#include <vector>

#include <unistd.h>
#include <sys/wait.h>

#include "larch/dag_loader.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"

#include "larch/usher_glue.hpp"
#include <tbb/task_scheduler_init.h>

template <typename Node1, typename Node2>
void compareDAG(Node1 dag1, Node2 dag2) {
  if (dag1.GetCladesCount() != dag2.GetCladesCount()) {
    fprintf(stderr, "Node %zu have %zu clades,but Node %zu have %zu\n",
            dag1.GetId().value, dag1.GetCladesCount(), dag2.GetId().value,
            dag2.GetCladesCount());
  }
  for (size_t child_idx = 0; child_idx < dag1.GetCladesCount(); child_idx++) {
    auto edge1 = dag1.GetClade(CladeIdx{child_idx}).at(0);
    auto edge2 = dag2.GetClade(CladeIdx{child_idx}).at(0);
    if (edge1.GetEdgeMutations() != edge2.GetEdgeMutations()) {
      fprintf(stderr, "edge %zu and edge %zu  have mismatch mutation\n",
              edge1.GetId().value, edge2.GetId().value);
    }
    compareDAG(edge1.GetChild(), edge2.GetChild());
  }
}

template <typename DAG>
void check_MAT_MADAG_Eq(MAT::Tree& tree, DAG init) {
  auto converted_dag = AddMATConversion(MADAGStorage{{}});
  converted_dag.View().BuildFromMAT(tree, init.GetReferenceSequence());
  compareDAG(converted_dag.View().GetRoot(), init.GetRoot());
}

template <typename DAG, typename RadiusCallback, typename ReassignCallback>
auto optimize_dag_direct(DAG dag, Move_Found_Callback& callback,
                         RadiusCallback&& radius_callback,
                         ReassignCallback&& reassign_callback) {
  static_assert(DAG::template contains_element_feature<Component::Node, MATConversion>);
  auto& tree = dag.GetMutableMAT();

  // TODO only in debug mode
  Mutation_Annotated_Tree::save_mutation_annotated_tree(tree, "before_optimize.pb");
  check_MAT_MADAG_Eq(tree, dag);

  Original_State_t origin_states;
  check_samples(tree.root, origin_states, &tree);
  reassign_states(tree, origin_states);
  reassign_callback.OnReassignedStates(tree);
  radius_callback(tree);

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
  // NOLINTNEXTLINE
  std::chrono::steady_clock::time_point end_time = start_time + std::chrono::hours(8);
  size_t ddepth = tree.get_max_level() * 2;
  std::cout << "maximum radius is " << std::to_string(ddepth) << "\n";
  size_t rad_exp = 1;
  for (; static_cast<size_t>(1) << rad_exp <= ddepth; rad_exp++) {
    auto all_nodes = tree.depth_first_expansion();
    std::cout << "current radius is " << std::to_string(1 << rad_exp) << "\n";

    optimize_inner_loop(all_nodes,     // nodes to search
                        tree,          // tree
                        1 << rad_exp,  // radius
                        callback,
                        true,                  // allow drift
                        true,                  // search all directions
                        5,                     // NOLINT // minutes between save
                        true,                  // do not write intermediate files
                        end_time,              // search end time
                        start_time,            // start time
                        false,                 // log moves
                        1,                     // current iteration
                        "intermediate",        // intermediate template
                        "intermediate_base",   // intermediate base name
                        "intermediate_newick"  // intermediate newick name
    );
    radius_callback(tree);
  }
  Mutation_Annotated_Tree::save_mutation_annotated_tree(tree, "after_optimize.pb");
  auto result = std::make_pair(AddMATConversion(MADAGStorage{{}}), std::move(tree));
  result.first.View().BuildFromMAT(result.second, dag.GetReferenceSequence());
  // TODO tree.delete_nodes();
  result.first.View().RecomputeCompactGenomes(true);
  return result;
}

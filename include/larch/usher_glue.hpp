#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-user-provided-copy"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wconversion"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "src/matOptimize/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#pragma GCC diagnostic pop

#include "larch/madag/mutation_annotated_dag.hpp"

namespace Mutation_Annotated_Tree {
class Tree;
class Node;
}  // namespace Mutation_Annotated_Tree
namespace MAT = Mutation_Annotated_Tree;

inline void fill_static_reference_sequence(std::string_view dag_ref);

template <typename DAG>
void check_edge_mutations(DAG dag);

template <typename DAG>
MAT::Tree mat_from_dag(DAG dag);

template <typename DAG>
MADAGStorage optimize_dag_direct(DAG dag, Move_Found_Callback& callback);

#include "larch/impl/usher_glue_impl.hpp"
#include "larch/impl/produce_mat_impl.hpp"
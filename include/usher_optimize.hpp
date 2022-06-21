#pragma once

#include <tbb/task_scheduler_init.h>
#include <mpi.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wparentheses"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wpedantic"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/Fitch_Sankoff.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#pragma GCC diagnostic pop

void InitUsherMPI(int argc, char **argv);

void UsherOptimize(Mutation_Annotated_Tree::Tree &t);

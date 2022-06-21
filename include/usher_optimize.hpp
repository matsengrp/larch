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

namespace MAT = Mutation_Annotated_Tree;

extern std::atomic_bool interrupted;
extern bool use_bound;
extern int process_count;
extern int this_rank;
extern uint32_t num_threads;

#define DRIFT_MASK 0x80000000
#define ALL_DIR_MASK 0x40000000

int count_back_mutation(const MAT::Tree &tree);
void get_pos_samples_old_tree(MAT::Tree &tree, std::vector<mutated_t> &output);
void min_back_reassign_state_local(MAT::Tree &tree,
                                   const std::vector<mutated_t> &mutations);
void MPI_min_back_reassign_states(MAT::Tree &tree,
                                  const std::vector<mutated_t> &mutations,
                                  int start_position);

void UsherOptimize(Mutation_Annotated_Tree::Tree &t);

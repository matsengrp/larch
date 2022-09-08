#include <cstddef>
#include <iostream>
#include <charconv>

#include <unistd.h>
#include <sys/wait.h>

#include "arguments.hpp"
#include "dag_loader.hpp"
#include "edge.hpp"
#include "subtree_weight.hpp"
#include "parsimony_score.hpp"
#include "merge.hpp"

#include "../deps/usher/src/matOptimize/mutation_annotated_tree.hpp"
#include "../deps/usher/src/matOptimize/check_samples.hpp"
#include "../deps/usher/src/matOptimize/tree_rearrangement_internal.hpp"
#include <tbb/task_scheduler_init.h>

namespace MAT = Mutation_Annotated_Tree;
std::atomic_bool interrupted(false);
int process_count=1;
int this_rank=0;
uint32_t num_threads=tbb::task_scheduler_init::default_num_threads();
FILE* movalbe_src_log;
bool changing_radius=false;
static int32_t EncodeBase(char base) {
  switch (base) {
    case 'A':
      return 1;
    case 'C':
      return 2;
    case 'G':
      return 4;
    case 'T':
      return 8;
    default:
      Fail("Invalid base");
  };
}
static void mat_from_dag_helper(Node dag_node, MAT::Node* mat_par_node,const std::vector<EdgeMutations>& edge_mutations,MAT::Tree& new_tree ){
    mat_par_node->children.reserve(dag_node.GetCladesCount());
    for(auto dag_node_edge_child:dag_node.GetClades()){
        auto edge=dag_node_edge_child[0];
        const auto& mutations=edge_mutations[edge.GetId().value];
        auto node=new MAT::Node(edge.GetChildId().value);
        new_tree.register_node_serial(node);
        node->mutations.reserve(mutations.size());
        for (const auto& mut: mutations){
            MAT::Mutation mat_mut("ref",mut.first.value,EncodeBase(mut.second.second),EncodeBase(mut.second.first),EncodeBase(mut.second.second));
            node->mutations.push_back(mat_mut);
        }
        node->parent=mat_par_node;
        mat_par_node->children.push_back(node);
        mat_from_dag_helper(edge.GetChild(), node, edge_mutations, new_tree);
    }
}
bool use_bound=true;


MAT::Tree mat_from_dag(const MADAG& dag){
    MAT::Tree tree;
    const std::vector<EdgeMutations>& edge_mutations=dag.GetEdgeMutations();
    auto root_node=dag.GetDAG().GetRoot();
    MAT::Node* mat_root_node=new MAT::Node (root_node.GetId().value);
    tree.root=mat_root_node;
    tree.register_node_serial(mat_root_node);
    mat_from_dag_helper(root_node, mat_root_node, edge_mutations, tree);
    return tree;
}

void check_mat_madag_conversion(const MADAG& dag){
    auto dag_ref=dag.GetReferenceSequence();
    MAT::Mutation::refs.resize(dag_ref.size());
    for (size_t ref_idx=0; ref_idx<dag_ref.size(); ref_idx++) {
        MAT::Mutation::refs[ref_idx]=EncodeBase(dag_ref[ref_idx]);
    }
    auto tree=mat_from_dag(dag);
    Original_State_t origin_states;
    check_samples(tree.root, origin_states, &tree);
    reassign_states(tree, origin_states);
    auto all_nodes=tree.depth_first_expansion();
    optimize_inner_loop(all_nodes, tree, 2);
    tree.delete_nodes();
}

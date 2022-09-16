#include <cstddef>
#include <cstdio>
#include <iostream>
#include <charconv>
#include <vector>

#include <unistd.h>
#include <sys/wait.h>

#include "arguments.hpp"
#include "dag.hpp"
#include "dag_loader.hpp"
#include "edge.hpp"
#include "edge_mutations.hpp"
#include "mutation_annotated_dag.hpp"
#include "node.hpp"
#include "range/v3/view/transform.hpp"
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
void  build_madag_from_mat_helper(MAT::Node* par_node,size_t par_node_idx, size_t& curr_idx,DAG& dag,std::vector<EdgeMutations>& edge_muts,size_t& edge_idx){
  for (size_t clade_idx=0; clade_idx<par_node->children.size(); clade_idx++) {
    curr_idx++;
    dag.AddEdge(EdgeId{edge_idx} ,NodeId{par_node_idx}, NodeId{curr_idx}, CladeIdx{clade_idx});
    const auto child_node=par_node->children[clade_idx];
    auto mut_view=child_node->mutations| ranges::views::transform(
                 [](const MAT::Mutation& mut) -> std::pair<MutationPosition, std::pair<char, char>> {
                   static const char decode[] = {'A', 'C', 'G', 'T'};
                   return {{static_cast<size_t>(mut.get_position())},
                           {decode[one_hot_to_two_bit(mut.get_par_one_hot())],
                            decode[one_hot_to_two_bit(mut.get_mut_one_hot())]}};
                 });
    edge_muts[edge_idx]=EdgeMutations(mut_view);
    edge_idx++;
    build_madag_from_mat_helper(child_node, curr_idx, curr_idx, dag, edge_muts, edge_idx);
  }
} 
MADAG build_madag_from_mat(const MAT::Tree& tree){
  MADAG result;
  auto & edge_muts=result.GetEdgeMutations();
  edge_muts.resize(tree.get_size_upper());
  size_t node_count=0;
  size_t edge_count=0;
  build_madag_from_mat_helper(tree.root, 0, node_count, result.GetDAG(), edge_muts, edge_count);
  result.GetDAG().InitializeNodes(node_count+1);
  result.GetDAG().BuildConnections();
  edge_muts.resize(edge_count+1);
  return result;
}
void compareDAG(const Node dag1,const Node dag2,const std::vector<EdgeMutations>& edge_mutations1,const std::vector<EdgeMutations>& edge_mutations2 ){
  if(dag1.GetCladesCount()!=dag2.GetCladesCount()){
    fprintf(stderr, "Node %zu have %zu clades,but Node %zu have %zu\n", dag1.GetId().value,dag1.GetCladesCount(),dag2.GetId().value,dag2.GetCladesCount());
  }
  for(size_t child_idx=0;child_idx<dag1.GetCladesCount();child_idx++ ){
    auto edge1=dag1.GetClade(CladeIdx{child_idx})[0];
    auto edge2=dag2.GetClade(CladeIdx{child_idx})[0];
    if(edge_mutations1[edge1.GetId().value]!=edge_mutations2[edge1.GetId().value]){
      fprintf(stderr, "edge %zu and edge %zu  have mismatch mutation\n",edge1.GetId().value,edge2.GetId().value);
    }
    compareDAG(edge1.GetChild(), edge2.GetChild(),edge_mutations1,edge_mutations2);
  }
}
void check_MAT_MADAG_Eq(const MAT::Tree& tree,const MADAG& init){
  MADAG converted_dag =build_madag_from_mat(tree);
  compareDAG(converted_dag.GetDAG().GetRoot(), init.GetDAG().GetRoot(),converted_dag.GetEdgeMutations(),converted_dag.GetEdgeMutations());
}
MADAG optimize_dag_direct(const MADAG& dag){
    auto dag_ref=dag.GetReferenceSequence();
    MAT::Mutation::refs.resize(dag_ref.size()+1);
    for (size_t ref_idx=0; ref_idx<dag_ref.size(); ref_idx++) {
        MAT::Mutation::refs[ref_idx+1]=EncodeBase(dag_ref[ref_idx]);
    }
    auto tree=mat_from_dag(dag);
    check_MAT_MADAG_Eq(tree, dag);
    Original_State_t origin_states;
    check_samples(tree.root, origin_states, &tree);
    reassign_states(tree, origin_states);
    auto all_nodes=tree.depth_first_expansion();

    std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point end_time = start_time + std::chrono::hours(8);
    int ddepth = tree.get_max_level() * 2;
    std::cout << "maximum radius is " << std::to_string(ddepth) << "\n";
    for (size_t rad_exp = 1; (1 << rad_exp) <= ddepth; rad_exp++) {
        std::cout << "current radius is " << std::to_string(1 << rad_exp) << "\n";

        optimize_inner_loop(all_nodes, // nodes to search
                tree, // tree
                1 << rad_exp, // radius
                true, // allow drift
                true, // search all directions
                5, // minutes between save
                true, // do not write intermediate files
                end_time, // search end time
                start_time, // start time
                false, // log moves
                1, // current iteration
                "intermediate", // intermediate template
                "intermediate_base", // intermediate base name
                "intermediate_newick" // intermediate newick name
                );
    }
    MADAG result=build_madag_from_mat(tree);
    tree.delete_nodes();
    return result;
}

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <string_view>
#include <vector>

#include <unistd.h>
#include <sys/wait.h>

#include "larch/dag/dag.hpp"
#include "larch/dag_loader.hpp"
#include "larch/dag/edge.hpp"
#include "larch/merge/edge_mutations.hpp"
#include "larch/mutation_annotated_dag.hpp"
#include "larch/dag/node.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/merge/merge.hpp"

#include "larch/usher_glue.hpp"
#include <tbb/task_scheduler_init.h>

namespace MAT = Mutation_Annotated_Tree;
std::atomic_bool interrupted(false);
int process_count = 1;
int this_rank = 0;
uint32_t num_threads =
    static_cast<uint32_t>(tbb::task_scheduler_init::default_num_threads());
FILE* movalbe_src_log;
bool changing_radius = false;
bool use_bound = true;

static uint8_t EncodeBase(char base) {
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
static void mat_from_dag_helper(Node dag_node, MAT::Node* mat_par_node, size_t& node_id,
                                const std::vector<EdgeMutations>& edge_mutations,
                                MAT::Tree& new_tree) {
  mat_par_node->children.reserve(dag_node.GetCladesCount());
  for (auto clade : dag_node.GetClades()) {
    Assert(clade.size() == 1);
    Edge edge = *clade.begin();
    const auto& mutations = edge_mutations.at(edge.GetId().value);
    MAT::Node* node = new MAT::Node(node_id++);
    new_tree.register_node_serial(node);
    node->mutations.reserve(mutations.size());
    for (auto [pos, muts] : mutations) {
      Assert(pos.value != NoId);
      MAT::Mutation mat_mut("ref", static_cast<int>(pos.value), EncodeBase(muts.second),
                            EncodeBase(muts.first), EncodeBase(muts.second));
      node->mutations.push_back(mat_mut);
    }
    node->parent = mat_par_node;
    mat_par_node->children.push_back(node);
    mat_from_dag_helper(edge.GetChild(), node, node_id, edge_mutations, new_tree);
  }
}

MAT::Tree mat_from_dag(const MADAG& dag) {
  dag.AssertUA();
  MAT::Tree tree;
  size_t node_id = 0;
  const std::vector<EdgeMutations>& edge_mutations = dag.GetEdgeMutations();
  Node root_node = dag.GetDAG().GetRoot().GetFirstChild().GetChild();
  MAT::Node* mat_root_node = new MAT::Node(node_id++);

  const auto& tree_root_mutations =
      edge_mutations.at(dag.GetDAG().GetRoot().GetFirstChild().GetId().value);
  mat_root_node->mutations.reserve(tree_root_mutations.size());
  for (auto [pos, muts] : tree_root_mutations) {
    Assert(pos.value != NoId);
    MAT::Mutation mat_mut("ref", static_cast<int>(pos.value), EncodeBase(muts.second),
                          EncodeBase(muts.first), EncodeBase(muts.second));
    mat_root_node->mutations.push_back(mat_mut);
  }

  tree.root = mat_root_node;
  tree.register_node_serial(mat_root_node);
  mat_from_dag_helper(root_node, mat_root_node, node_id, edge_mutations, tree);

  return tree;
}

auto mutations_view(MAT::Node* node) {
  return node->mutations |
         ranges::views::transform(
             [](const MAT::Mutation& mut)
                 -> std::pair<MutationPosition, std::pair<char, char>> {
               static const char decode[] = {'A', 'C', 'G', 'T'};
               return {{static_cast<size_t>(mut.get_position())},
                       {decode[one_hot_to_two_bit(mut.get_par_one_hot())],
                        decode[one_hot_to_two_bit(mut.get_mut_one_hot())]}};
             });
}

void build_madag_from_mat_helper(MAT::Node* par_node, Node node, MADAG& dag) {
  for (size_t clade_idx = 0; clade_idx < par_node->children.size(); clade_idx++) {
    MAT::Node* mat_child = par_node->children[clade_idx];
    Node child_node = dag.AppendNode();
    Edge child_edge = dag.AppendEdge(node, child_node, CladeIdx{clade_idx});
    dag.ResizeEdgeMutations(dag.GetDAG().GetEdgesCount());
    dag.SetEdgeMutations(child_edge, EdgeMutations{mutations_view(mat_child)});
    build_madag_from_mat_helper(mat_child, child_node, dag);
  }
}

MADAG build_madag_from_mat(const MAT::Tree& tree, std::string_view reference_sequence) {
  MADAG result{reference_sequence};
  Node root_node = result.AppendNode();
  build_madag_from_mat_helper(tree.root, root_node, result);
  result.BuildConnections();
  result.AddUA(EdgeMutations{mutations_view(tree.root)});
  return result;
}

void compareDAG(const Node dag1, const Node dag2,
                const std::vector<EdgeMutations>& edge_mutations1,
                const std::vector<EdgeMutations>& edge_mutations2) {
  if (dag1.GetCladesCount() != dag2.GetCladesCount()) {
    fprintf(stderr, "Node %zu have %zu clades,but Node %zu have %zu\n",
            dag1.GetId().value, dag1.GetCladesCount(), dag2.GetId().value,
            dag2.GetCladesCount());
  }
  for (size_t child_idx = 0; child_idx < dag1.GetCladesCount(); child_idx++) {
    auto edge1 = dag1.GetClade(CladeIdx{child_idx})[0];
    auto edge2 = dag2.GetClade(CladeIdx{child_idx})[0];
    if (edge_mutations1[edge1.GetId().value] != edge_mutations2[edge2.GetId().value]) {
      fprintf(stderr, "edge %zu and edge %zu  have mismatch mutation\n",
              edge1.GetId().value, edge2.GetId().value);
    }
    compareDAG(edge1.GetChild(), edge2.GetChild(), edge_mutations1, edge_mutations2);
  }
}
void check_MAT_MADAG_Eq(const MAT::Tree& tree, const MADAG& init) {
  MADAG converted_dag = build_madag_from_mat(tree, init.GetReferenceSequence());
  compareDAG(converted_dag.GetDAG().GetRoot(), init.GetDAG().GetRoot(),
             converted_dag.GetEdgeMutations(), init.GetEdgeMutations());
}

void fill_static_reference_sequence(std::string_view dag_ref) {
  MAT::Mutation::refs.resize(dag_ref.size() + 1);
  for (size_t ref_idx = 0; ref_idx < dag_ref.size(); ref_idx++) {
    MAT::Mutation::refs[ref_idx + 1] = EncodeBase(dag_ref[ref_idx]);
  }
}

MADAG optimize_dag_direct(const MADAG& dag, Move_Found_Callback& callback) {
  auto dag_ref = dag.GetReferenceSequence();
  fill_static_reference_sequence(dag_ref);
  auto tree = mat_from_dag(dag);

  Mutation_Annotated_Tree::save_mutation_annotated_tree(tree, "before_optimize.pb");
  check_MAT_MADAG_Eq(tree, dag);
  Original_State_t origin_states;
  check_samples(tree.root, origin_states, &tree);
  reassign_states(tree, origin_states);

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point end_time = start_time + std::chrono::hours(8);
  size_t ddepth = tree.get_max_level() * 2;
  std::cout << "maximum radius is " << std::to_string(ddepth) << "\n";
  size_t rad_exp = 1;
  // for (; static_cast<size_t>(1) << rad_exp <= ddepth; rad_exp++) {
    auto all_nodes = tree.depth_first_expansion();
    std::cout << "current radius is " << std::to_string(1 << rad_exp) << "\n";

    optimize_inner_loop(all_nodes,     // nodes to search
                        tree,          // tree
                        1 << rad_exp,  // radius
                        callback,
                        true,                  // allow drift
                        true,                  // search all directions
                        5,                     // minutes between save
                        true,                  // do not write intermediate files
                        end_time,              // search end time
                        start_time,            // start time
                        false,                 // log moves
                        1,                     // current iteration
                        "intermediate",        // intermediate template
                        "intermediate_base",   // intermediate base name
                        "intermediate_newick"  // intermediate newick name
    );
  // }
  Mutation_Annotated_Tree::save_mutation_annotated_tree(tree, "after_optimize.pb");
  MADAG result = build_madag_from_mat(tree, dag.GetReferenceSequence());
  tree.delete_nodes();
  return result;
}

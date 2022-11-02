#include "larch/dag_loader.hpp"
#include "larch/dag/edge.hpp"
#include "larch/merge/edge_mutations.hpp"
#include "larch/mutation_annotated_dag.hpp"
#include "larch/dag/node.hpp"
#include <iostream>
#include <ostream>
#include <unordered_map>
#include <vector>
struct ParNuc_Info {
  NodeId par_node_id;
  NodeId child_node_id;
  EdgeId edge_id;
  char nuc;
};
typedef std::unordered_map<size_t, ParNuc_Info> all_mutated_t;
void check_edge_mutations_helper(const Node dag_node, const all_mutated_t& all_mutated,
                                 const std::string_view ref_seq,
                                 const std::vector<EdgeMutations>& edge_mutations) {
  for (auto child : dag_node.GetChildren()) {
    auto edge_id = child.GetId();
    all_mutated_t this_mutated(all_mutated);
    for (const auto& edge_mut : edge_mutations[edge_id.value]) {
      auto ins_result = this_mutated.emplace(
          edge_mut.first.value, ParNuc_Info{dag_node.GetId(), child.GetChildId(),
                                            edge_id, edge_mut.second.second});
      auto actual_par = edge_mut.second.first;
      if (ins_result.second) {
        auto expected = ref_seq[edge_mut.first.value - 1];
        if (expected != actual_par) {
          std::cout << "On edge " << edge_id.value << " from " << dag_node.GetId().value
                    << " to " << child.GetChildId().value << " at position "
                    << edge_mut.first.value << " expected reference " << expected
                    << " but got " << actual_par << "\n";
        }
      } else {
        auto expected = ins_result.first->second.nuc;
        if (expected != actual_par) {
          std::cout << "On edge " << edge_id.value << " from " << dag_node.GetId().value
                    << " to " << child.GetChildId().value << " at position "
                    << edge_mut.first.value << " expected " << expected
                    << " from mutation on edge "
                    << ins_result.first->second.edge_id.value << "from "
                    << ins_result.first->second.par_node_id.value << " to "
                    << ins_result.first->second.child_node_id.value << " but got "
                    << actual_par << "\n";
        }
        ins_result.first->second.child_node_id = child.GetChildId();
        ins_result.first->second.par_node_id = dag_node.GetId();
        ins_result.first->second.edge_id = child.GetId();
        ins_result.first->second.nuc = edge_mut.second.second;
      }
    }
    check_edge_mutations_helper(child.GetChild(), this_mutated, ref_seq,
                                edge_mutations);
  }
}

void check_edge_mutations(const MADAG& madag) {
  std::cout << "start_check" << std::endl;
  const auto ref_req = madag.GetReferenceSequence();
  all_mutated_t init;
  check_edge_mutations_helper(madag.GetDAG().GetRoot(), init,
                              madag.GetReferenceSequence(), madag.GetEdgeMutations());
  std::cout << "end_check" << std::endl;
}
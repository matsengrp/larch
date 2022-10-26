#include "test_common.hpp"

#include "dag_loader.hpp"
#include "merge.hpp"
#include "subtree_weight.hpp"
#include "parsimony_score.hpp"

#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"

namespace MAT = Mutation_Annotated_Tree;

void check_edge_mutations(const MADAG&);
MADAG optimize_dag_direct(const MADAG&, Move_Found_Callback&);

struct Test_Move_Found_Callback : public Move_Found_Callback {
  bool operator()(Profitable_Moves& move, int best_score_change,
                  std::vector<Node_With_Major_Allele_Set_Change>&) override {
    return move.score_change < best_score_change;
  }
};

static void test_matOptimize(std::string_view input_dag_path,
                             std::string_view refseq_path, size_t count) {
  std::string reference_sequence = LoadReferenceSequence(refseq_path);
  MADAG input_dag = LoadTreeFromProtobuf(input_dag_path, reference_sequence);
  Merge merge{input_dag.GetReferenceSequence()};
  merge.AddDAGs({input_dag});
  std::vector<MADAG> optimized_dags;

  for (size_t i = 0; i < count; ++i) {
    merge.ComputeResultEdgeMutations();
    SubtreeWeight<ParsimonyScore> weight{merge.GetResult()};
    auto [sample, dag_ids] = weight.SampleTree({});
    check_edge_mutations(sample);
    MADAG result;
    Test_Move_Found_Callback callback;
    result = optimize_dag_direct(sample, callback);
    optimized_dags.push_back(std::move(result));
    merge.AddDAGs({optimized_dags.back()});
  }
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] {
                test_matOptimize("data/startmat/startmat_no_ancestral.pb.gz",
                                 "data/startmat/refseq.txt.gz", 3);
              },
              "matOptimize: tree startmat"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] {
                test_matOptimize("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
                                 "data/20D_from_fasta/refseq.txt.gz", 3);
              },
              "matOptimize: tree 20D_from_fasta"});

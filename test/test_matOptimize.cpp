#include <tuple>
#include "test_common.hpp"

#include "larch/dag_loader.hpp"
#include "larch/merge/merge.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"

#include "larch/usher_glue.hpp"

namespace MAT = Mutation_Annotated_Tree;

void check_edge_mutations(MADAG);
MADAGStorage optimize_dag_direct(MADAG, Move_Found_Callback&);

struct Test_Move_Found_Callback : public Move_Found_Callback {
  bool operator()(Profitable_Moves& move, int best_score_change,
                  std::vector<Node_With_Major_Allele_Set_Change>&) override {
    return move.score_change < best_score_change;
  }
};

static void test_matOptimize(std::string_view input_dag_path,
                             std::string_view refseq_path, size_t count) {
  std::string reference_sequence = LoadReferenceSequence(refseq_path);
  MADAGStorage input_dag_storage =
      LoadTreeFromProtobuf(input_dag_path, reference_sequence);
  MADAG input_dag = input_dag_storage.View();
  Merge merge{input_dag.GetReferenceSequence()};
  merge.AddDAGs({input_dag});
  std::vector<MADAGStorage> optimized_dags;

  for (size_t i = 0; i < count; ++i) {
    merge.ComputeResultEdgeMutations();
    SubtreeWeight<ParsimonyScore> weight{merge.GetResult()};
    auto [sample, dag_ids] = weight.SampleTree({});
    std::ignore = dag_ids;
    check_edge_mutations(sample.View());
    Test_Move_Found_Callback callback;
    optimized_dags.push_back(optimize_dag_direct(sample.View(), callback));
    merge.AddDAGs({optimized_dags.back().View()});
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

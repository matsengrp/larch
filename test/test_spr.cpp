#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/merge/merge.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/spr/spr_view.hpp"

template <typename DAG>
struct Test_Move_Found_Callback : public Move_Found_Callback {
  Test_Move_Found_Callback(DAG sample_dag) : sample_dag_{sample_dag} {}

  bool operator()(Profitable_Moves& move, int best_score_change,
                  std::vector<Node_With_Major_Allele_Set_Change>&
                      nodes_with_major_allele_set_change) override {
    auto spr_storage = SPRStorage(sample_dag_);
    auto spr = spr_storage.View();

    Assert(sample_mat_);
    spr.InitHypotheticalTree(sample_dag_, *sample_mat_, move,
                             nodes_with_major_allele_set_change);

    // std::ignore = spr.GetRoot().ComputeNewCompactGenome();

    return move.score_change < best_score_change;
  }

  void operator()(const MAT::Tree& tree) { sample_mat_ = std::addressof(tree); }

  DAG sample_dag_;
  const MAT::Tree* sample_mat_ = nullptr;
};

static void test_spr(std::string_view input_dag_path, std::string_view refseq_path,
                     size_t count) {
  std::string reference_sequence = LoadReferenceSequence(refseq_path);
  MADAGStorage input_dag_storage =
      LoadTreeFromProtobuf(input_dag_path, reference_sequence);
  input_dag_storage.View().RecomputeCompactGenomes();
  MADAG input_dag = input_dag_storage.View();
  Merge<MADAG> merge{input_dag.GetReferenceSequence()};
  merge.AddDAGs({input_dag});
  std::vector<decltype(AddMATConversion(MADAGStorage{}))> optimized_dags;

  for (size_t i = 0; i < count; ++i) {
    merge.ComputeResultEdgeMutations();
    SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};

    auto chosen_node = weight.GetDAG().GetRoot();
    auto sample = weight.SampleTree({}, chosen_node);
    std::cout << "Sample nodes count: " << sample.GetNodesCount() << "\n";
    check_edge_mutations(sample.View());
    Test_Move_Found_Callback callback{sample.View()};
    optimized_dags.push_back(optimize_dag_direct(sample.View(), callback, callback));
    optimized_dags.back().View().RecomputeCompactGenomes();
    merge.AddDAG(optimized_dags.back().View(), chosen_node);
  }
}

[[maybe_unused]] static const auto test_added1 =
    add_test({[] {
                test_spr("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
                         "data/20D_from_fasta/refseq.txt.gz", 3);
              },
              "SPR: tree 20D_from_fasta"});

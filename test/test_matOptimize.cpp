#include <tuple>
#include "test_common.hpp"

#include "larch/dag_loader.hpp"
#include "larch/merge/merge.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"

#include "larch/usher_glue.hpp"

struct Test_Move_Found_Callback : public Move_Found_Callback {
  bool operator()(Profitable_Moves& move, int best_score_change,
                  std::vector<Node_With_Major_Allele_Set_Change>&) override {
    return move.score_change < best_score_change;
  }
};

[[maybe_unused]] static auto choose_root = [](const auto& subtree_weight) {
  return subtree_weight.GetDAG().GetRoot();
};

static auto choose_random = [](const auto& weight) {
  std::random_device random_device;
  std::mt19937 random_generator(random_device());

  NodeId node_id;
  do {
    Assert(weight.GetDAG().GetNodesCount() > 0);
    node_id = {std::uniform_int_distribution<size_t>{
        0, weight.GetDAG().GetNodesCount() - 1}(random_generator)};
    auto node = weight.GetDAG().Get(node_id);
    if (node.IsLeaf() or node.GetCompactGenome().empty()) {
      continue;
    }
    break;
  } while (true);
  return weight.GetDAG().Get(node_id);
};

template <typename ChooseNode>
static void test_matOptimize(std::string_view input_dag_path,
                             std::string_view refseq_path, size_t count,
                             ChooseNode& choose_node) {
  std::string reference_sequence = LoadReferenceSequence(refseq_path);
  MADAGStorage input_dag_storage =
      LoadTreeFromProtobuf(input_dag_path, reference_sequence);
  input_dag_storage.View().RecomputeCompactGenomes();
  MADAG input_dag = input_dag_storage.View();
  Merge<MADAG> merge{input_dag.GetReferenceSequence()};
  merge.AddDAGs({input_dag});
  std::vector<MADAGStorage> optimized_dags;

  for (size_t i = 0; i < count; ++i) {
    merge.ComputeResultEdgeMutations();
    SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};

    auto chosen_node = choose_node(weight);
    auto [sample, dag_ids] = weight.SampleTree({}, chosen_node);
    std::cout << "Sample nodes count: " << sample.GetNodesCount() << "\n";
    std::ignore = dag_ids;
    check_edge_mutations(sample.View());
    Test_Move_Found_Callback callback;
    optimized_dags.push_back(optimize_dag_direct(sample.View(), callback));
    optimized_dags.back().View().RecomputeCompactGenomes();
    merge.AddDAG(optimized_dags.back().View(), chosen_node);
  }
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] {
                test_matOptimize("data/startmat/startmat_no_ancestral.pb.gz",
                                 "data/startmat/refseq.txt.gz", 100, choose_random);
              },
              "matOptimize: tree startmat"});

[[maybe_unused]] static const auto test_added1 = add_test(
    {[] {
       test_matOptimize("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
                        "data/20D_from_fasta/refseq.txt.gz", 100, choose_random);
     },
     "matOptimize: tree 20D_from_fasta"});

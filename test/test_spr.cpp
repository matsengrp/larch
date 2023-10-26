#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/spr/batching_callback.hpp"
#include "sample_dag.hpp"
#include "benchmark.hpp"

#include <tbb/global_control.h>

struct Empty_Callback : public Move_Found_Callback {
  bool operator()(Profitable_Moves& move, int best_score_change,
                  [[maybe_unused]] std::vector<Node_With_Major_Allele_Set_Change>&
                      nodes_with_major_allele_set_change) override {
    return move.score_change < best_score_change;
  }
  void operator()(MAT::Tree&) {}
  void OnReassignedStates(MAT::Tree&) {}
  auto GetMATNodeToCGMap() { return std::map<MAT::Node*, CompactGenome>{}; }
};

template <typename SampleDAG>
struct Test_Move_Found_Callback
    : public BatchingCallback<Test_Move_Found_Callback<SampleDAG>, SampleDAG> {
  explicit Test_Move_Found_Callback(Merge& merge, SampleDAG sample_dag)
      : BatchingCallback<Test_Move_Found_Callback<SampleDAG>, SampleDAG>{merge,
                                                                         sample_dag} {};

  template <typename SPRView, typename FragmentType>
  std::pair<bool, bool> OnMove(SPRView spr, const FragmentType& fragment,
                               Profitable_Moves& move, int best_score_change,
                               std::vector<Node_With_Major_Allele_Set_Change>&
                                   nodes_with_major_allele_set_change) {
    std::ignore = spr;
    std::ignore = fragment;
    std::ignore = nodes_with_major_allele_set_change;
    std::ignore = best_score_change;
    return {move.score_change <= 0, move.score_change <= 0};
  }

  void OnRadius() {}
};

[[maybe_unused]] static MADAGStorage<> Load(std::string_view input_dag_path,
                                            std::string_view refseq_path) {
  std::string reference_sequence = LoadReferenceSequence(refseq_path);
  MADAGStorage<> input_dag_storage =
      LoadTreeFromProtobuf(input_dag_path, reference_sequence);
  input_dag_storage.View().RecomputeCompactGenomes(true);
  return input_dag_storage;
}

[[maybe_unused]] static MADAGStorage<> Load(std::string_view input_dag_path) {
  MADAGStorage<> input_dag_storage = LoadDAGFromProtobuf(input_dag_path);
  input_dag_storage.View().RecomputeCompactGenomes(true);
  input_dag_storage.View().SampleIdsFromCG();
  return input_dag_storage;
}

static void test_spr(const MADAGStorage<>& input_dag_storage, size_t count) {
  tbb::global_control c(tbb::global_control::max_allowed_parallelism, 1);
  MADAG input_dag = input_dag_storage.View();
  Merge merge{input_dag.GetReferenceSequence()};
  merge.AddDAGs(std::vector{input_dag});
  std::vector<
      std::pair<decltype(AddMATConversion(MADAGStorage<>::EmptyDefault())), MAT::Tree>>
      optimized_dags;

  for (size_t i = 0; i < count; ++i) {
    merge.ComputeResultEdgeMutations();
    SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};

    auto chosen_node = weight.GetDAG().GetRoot();
    auto sample = AddMATConversion(weight.SampleTree({}, chosen_node));
    MAT::Tree mat;
    sample.View().GetRoot().Validate(true);
    sample.View().BuildMAT(mat);
    sample.View().GetRoot().Validate(true);
    check_edge_mutations(sample.View().Const());
    Test_Move_Found_Callback callback{merge, sample.View()};
    // Empty_Callback callback;
    optimized_dags.push_back(
        optimize_dag_direct(sample.View(), callback, callback, callback));
    optimized_dags.back().first.View().RecomputeCompactGenomes();
    merge.AddDAGs(std::vector{optimized_dags.back().first.View()},
                  optimized_dags.back().first.View().GetRoot());
  }
}

template <typename SampleDAG>
struct Single_Move_Callback_With_Hypothetical_Tree : public Move_Found_Callback {
  Single_Move_Callback_With_Hypothetical_Tree(Merge& merge, SampleDAG sample)
      : merge_{merge}, sample_{sample}, approved_a_move_{false} {}

  bool operator()(Profitable_Moves& move, int /*best_score_change*/,
                  std::vector<Node_With_Major_Allele_Set_Change>&
                      nodes_with_major_allele_set_change) override {
    if (!approved_a_move_) {
      // apply move to merge object.

      auto storage = [this](std::string ref_seq) {
        std::unique_lock<std::mutex> lock{mutex_};
        Assert(sample_mat_ != nullptr);
        using Storage = MergeDAGStorage;
        auto mat_conv = AddMATConversion(Storage::EmptyDefault());
        mat_conv.View().BuildFromMAT(*sample_mat_, ref_seq);
        check_edge_mutations(mat_conv.View());
        mat_conv.View().RecomputeCompactGenomes(true);
        return AddSPRStorage(std::move(mat_conv));
      }(sample_.GetReferenceSequence());
      auto spr = storage.View();

      // ** create hypothetical tree
      spr.InitHypotheticalTree(move, nodes_with_major_allele_set_change);

      // ** build fragment
      auto spr_fragment = spr.MakeFragment();

      // set flag so we don't approve any more moves
      approved_a_move_ = true;

      // ** merge fragment into merge
      // TODO merge_.AddDAG(spr_fragment);

      // return true so we do apply this move.
      return true;

    } else {
      return !approved_a_move_;
    }
  }

  void operator()(MAT::Tree& tree) {
    std::unique_lock<std::mutex> lock{mutex_};
    sample_mat_ = std::addressof(tree);
  }

  void OnReassignedStates(const MAT::Tree& tree) {
    for (auto leaf_node : tree.get_leaves()) {
      auto new_cg =
          sample_.GetNodeFromMAT(sample_.GetMAT().get_node(leaf_node->node_id))
              .GetCompactGenome()
              .Copy();
      mat_node_to_cg_map_[leaf_node] = new_cg.Copy();
    }
  }

  auto GetMATNodeToCGMap() { return std::map<MAT::Node*, CompactGenome>{}; }

  Merge& merge_;
  SampleDAG sample_;
  std::map<MAT::Node*, CompactGenome> mat_node_to_cg_map_;
  std::mutex mutex_;
  MAT::Tree* sample_mat_ = nullptr;
  bool approved_a_move_;
};

[[maybe_unused]] static void test_optimizing_with_hypothetical_tree(
    const MADAGStorage<>& tree_shaped_dag) {
  // tbb::global_control c(tbb::global_control::max_allowed_parallelism, 1);
  // this test takes a tree and uses matOptimize to apply a single move.

  Merge dag_altered_in_callback{tree_shaped_dag.View().GetReferenceSequence()};
  dag_altered_in_callback.AddDAGs(std::vector{tree_shaped_dag.View()});
  dag_altered_in_callback.ComputeResultEdgeMutations();

  // sample tree
  SubtreeWeight<ParsimonyScore, MergeDAG> weight{dag_altered_in_callback.GetResult()};
  auto sample = AddMATConversion(weight.SampleTree({}));
  MAT::Tree mat;
  sample.View().BuildMAT(mat);
  check_edge_mutations(sample.View());

  // create a callback that only allows one move
  Single_Move_Callback_With_Hypothetical_Tree single_move_callback{
      dag_altered_in_callback, sample.View()};

  // optimize tree with matOptimize using a callback that only applies a single move
  auto [optimized_dag, optimized_mat] = optimize_dag_direct(
      sample.View(), single_move_callback, single_move_callback, single_move_callback);

  optimized_dag.View().RecomputeCompactGenomes(true);
  Merge two_tree_dag{tree_shaped_dag.View().GetReferenceSequence()};
  two_tree_dag.AddDAGs(std::vector{sample.View()});
  two_tree_dag.AddDAGs(std::vector{optimized_dag.View()});

  // check topologies of the two DAGs match in feature count
  Assert(two_tree_dag.GetResult().GetNodesCount() ==
         dag_altered_in_callback.GetResult().GetNodesCount());
  Assert(two_tree_dag.GetResult().GetEdgesCount() ==
         dag_altered_in_callback.GetResult().GetEdgesCount());
}

[[maybe_unused]] static void test_sample() {
  auto input_storage = make_sample_dag();
  auto dag = input_storage.View();
  auto spr_storage = AddSPRStorage(dag);
  auto spr = spr_storage.View();

  spr.GetRoot().Validate(true);

  for (auto node : spr.GetNodes()) {
    if (not node.IsOverlaid<CompactGenome>()) {
      node.SetOverlay<CompactGenome>();
    }
  }

  spr.RecomputeCompactGenomes(true);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] {
                auto input = Load("data/test_5_trees/tree_0.pb.gz");
                Benchmark bench;
                bench.start();
                test_spr(input, 3);
                bench.stop();
                std::cout << "Time: " << bench.durationMs() << " ms ";
              },
              "SPR: test_5_trees"});

// [[maybe_unused]] static const auto test_added1 =
//     add_test({[] {
//                 test_spr(Load("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
//                               "data/20D_from_fasta/refseq.txt.gz"),
//                          1);
//               },
//               "SPR: tree 20D_from_fasta"});

[[maybe_unused]] static const auto test_added2 =
    add_test({[] { test_spr(make_sample_dag(), 10); }, "SPR: sample"});

[[maybe_unused]] static const auto test_added3 =
    add_test({[] { test_sample(); }, "SPR: move"});

//[[maybe_unused]] static const auto test_added4 = add_test(
//    {[] {
//       test_spr(Load("data/seedtree/seedtree.pb.gz", "data/seedtree/refseq.txt.gz"),
//       4);
//     },
//     "SPR: seedtree"});

// [[maybe_unused]] static const auto test_added4 =
//     add_test({[] {
//                 test_optimizing_with_hypothetical_tree(
//                     Load("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
//                          "data/20D_from_fasta/refseq.txt.gz"));
//               },
//               "SPR: single"});

// [[maybe_unused]] static const auto test_added5 =
//     add_test({[] {
//                 test_optimizing_with_hypothetical_tree(
//                     Load("data/seedtree/seedtree.pb.gz",
//                     "data/seedtree/refseq.txt.gz")

//                 );
//               },
//               "SPR: seedtree single move"});

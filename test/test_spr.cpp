#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/spr/batching_callback.hpp"
#include "test_common_dag.hpp"
#include "larch/benchmark.hpp"

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

struct Test_Move_Found_Callback : public BatchingCallback<Test_Move_Found_Callback> {
  explicit Test_Move_Found_Callback(Merge& merge)
      : BatchingCallback<Test_Move_Found_Callback>{merge} {};

  template <typename SPRView, typename FragmentType>
  std::pair<bool, bool> OnMove(SPRView spr, const FragmentType& fragment,
                               Profitable_Moves& move, int best_score_change,
                               std::vector<Node_With_Major_Allele_Set_Change>&
                                   nodes_with_major_allele_set_change) {
    std::ignore = spr;
    std::ignore = fragment;
    std::ignore = nodes_with_major_allele_set_change;
    std::ignore = best_score_change;
    // for (auto hypothetical_node : fragment.GetNodes()) {
    //   std::cout << "in OnMove, node " << hypothetical_node << " has "
    //             << hypothetical_node.GetCladesCount() << " children\n"
    //             << std::flush;
    // }
    // for (auto node : spr.GetNodes()) {
    //   std::cout << "in OnMove, node " << node << " has " << node.GetCladesCount()
    //             << " children\n"
    //             << std::flush;
    // }
    // if (moves_count_.fetch_add(1) > 100) {
    //   print_peak_mem();
    //   moves_count_.store(0);
    // }
    // return {false, false};
    return {move.score_change <= 0, move.score_change <= 0};
  }

  void OnRadius() { std::cout << "OnRadius\n"; }
  std::atomic<size_t> moves_count_{0};
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
  input_dag_storage.View().SampleIdsFromCG(true);
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
    std::cout << "Computing edge mutations\n";
    merge.ComputeResultEdgeMutations();
    SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};

    auto chosen_node = weight.GetDAG().GetRoot();
    std::cout << "Sampling tree\n";
    auto sample = AddMATConversion(weight.SampleTree({}, chosen_node));
    MAT::Tree mat;
    sample.View().GetRoot().Validate(true);
    std::cout << "Building MAT\n";
    sample.View().BuildMAT(mat);
    sample.View().GetRoot().Validate(true);
    check_edge_mutations(sample.View().Const());
    Test_Move_Found_Callback callback{merge};
    std::cout << "Optimizing\n";
    // Empty_Callback callback;
    optimized_dags.push_back(
        optimize_dag_direct(sample.View(), callback, callback, callback));
    optimized_dags.back().first.View().RecomputeCompactGenomes();
    merge.AddDAGs(std::vector{optimized_dags.back().first.View()},
                  optimized_dags.back().first.View().GetRoot());
    mat.delete_nodes();
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
        TestAssert(sample_mat_ != nullptr);
        using Storage = MergeDAGStorage<>;
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
      auto node = sample_.GetNodeFromMAT(sample_.GetMAT().get_node(leaf_node->node_id));
      auto new_cg = node.GetCompactGenome().Copy(&node);
      mat_node_to_cg_map_[leaf_node] = std::move(new_cg);
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
  TestAssert(two_tree_dag.GetResult().GetNodesCount() ==
             dag_altered_in_callback.GetResult().GetNodesCount());
  TestAssert(two_tree_dag.GetResult().GetEdgesCount() ==
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

[[maybe_unused]] static const auto test_added1 =
    add_test({[] {
                auto input =
                    Load("data/seedtree/seedtree.pb.gz", "data/seedtree/refseq.txt.gz");
                Benchmark bench;
                bench.start();
                test_spr(input, 3);
                bench.stop();
                std::cout << "Time: " << bench.durationMs() << " ms ";
              },
              "SPR: seedtree",
              {"slow"}});

[[maybe_unused]] static const auto test_added2 =
    add_test({[] { test_spr(make_sample_dag(), 10); }, "SPR: sample"});

[[maybe_unused]] static const auto test_added3 =
    add_test({[] { test_sample(); }, "SPR: move"});

[[maybe_unused]] static const auto test_added4 =
    add_test({[] {
                Benchmark bench_load;
                std::cout << "Loading...\n";
                bench_load.start();
                auto input = Load("data/20B/20B_start_tree_no_ancestral.pb.gz",
                                  "data/20B/ref_seq_noancestral.txt.gz");
                bench_load.stop();
                std::cout << "Load time: " << bench_load.durationMs() << " ms\n";
                std::cout << "Node count: " << input.View().GetNodesCount() << "\n";
                Benchmark bench;
                bench.start();
                test_spr(input, 1);
                bench.stop();
                std::cout << "Time: " << bench.durationMs() << " ms ";
              },
              "SPR: 20B",
              {"slow"}});

[[maybe_unused]] static const auto test_added5 =
    add_test({[] {
                auto input = Load("data/ebov_dud17/output.pb.gz",
                                  "data/ebov_dud17/output.txt.gz");
                test_spr(input, 3);
              },
              "SPR: ebov_dud17",
              {"slow"}});

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

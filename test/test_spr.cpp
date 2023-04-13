#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/merge/merge.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/spr/spr_view.hpp"

#include <tbb/global_control.h>

template <typename DAG, typename MergeT>
struct Test_Move_Found_Callback : public Move_Found_Callback {
  Test_Move_Found_Callback(DAG sample_dag, MergeT& merge)
      : sample_dag_{sample_dag}, merge_{merge} {};

  using Storage =
      ExtendDAGStorage<DefaultDAGStorage,
                       Extend::Nodes<Deduplicate<CompactGenome>, SampleId>,
                       Extend::Edges<EdgeMutations>, Extend::DAG<ReferenceSequence>>;

  bool operator()(Profitable_Moves& move, int best_score_change,
                  [[maybe_unused]] std::vector<Node_With_Major_Allele_Set_Change>&
                      nodes_with_major_allele_set_change) override {
    Assert(move.src != nullptr);
    Assert(move.dst != nullptr);
    auto storage = [this](std::string ref_seq) {
      MAT::Tree* mat = sample_mat_.load();
      Assert(mat != nullptr);
      auto mat_conv = AddMATConversion(Storage{});
      mat_conv.View().BuildFromMAT(*mat, ref_seq);
      check_edge_mutations(mat_conv.View().Const());
      mat_conv.View().RecomputeCompactGenomes();
      return SPRStorage(std::move(mat_conv));
    }(sample_dag_.GetReferenceSequence());

    auto spr = storage.View();
    spr.GetRoot().Validate(true);
    if (spr.InitHypotheticalTree(move, nodes_with_major_allele_set_change)) {
      spr.GetRoot().Validate(true);

      std::cout << "After move: " << spr.GetMoveSource().GetId().value <<
        "->" << spr.GetMoveTarget().GetId().value << "\n";
      MADAGToDOT(spr, std::cout);
      
      auto fragment = spr.GetFragment();
      merge_.AddFragment(spr, fragment.first, fragment.second);
    } else {
      std::cout << "Skip move\n";
      return false;
    }
    return move.score_change < best_score_change;
  }

  void operator()(MAT::Tree& tree) { sample_mat_.store(std::addressof(tree)); }

  void OnReassignedStates(MAT::Tree& tree) {
    reassigned_states_storage_.View().BuildFromMAT(tree,
                                                   sample_dag_.GetReferenceSequence());
    check_edge_mutations(reassigned_states_storage_.View().Const());
    reassigned_states_storage_.View().RecomputeCompactGenomes();
    merge_.AddDAG(reassigned_states_storage_.View());
    merge_.ComputeResultEdgeMutations();
  }

  DAG sample_dag_;
  MergeT& merge_;
  decltype(AddMATConversion(Storage{})) reassigned_states_storage_ =
      AddMATConversion(Storage{});
  std::atomic<MAT::Tree*> sample_mat_ = nullptr;
};

static MADAGStorage Load(std::string_view input_dag_path,
                         std::string_view refseq_path) {
  std::string reference_sequence = LoadReferenceSequence(refseq_path);
  MADAGStorage input_dag_storage =
      LoadTreeFromProtobuf(input_dag_path, reference_sequence);
  input_dag_storage.View().RecomputeCompactGenomes();
  return input_dag_storage;
}

static void test_spr(const MADAGStorage& input_dag_storage, size_t count) {
  tbb::global_control c(tbb::global_control::max_allowed_parallelism, 1);
  MADAG input_dag = input_dag_storage.View();
  Merge<MADAG> merge{input_dag.GetReferenceSequence()};
  merge.AddDAG(input_dag);
  std::vector<std::pair<decltype(AddMATConversion(MADAGStorage{})), MAT::Tree>>
      optimized_dags;

  for (size_t i = 0; i < count; ++i) {
    merge.ComputeResultEdgeMutations();
    SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};

    auto chosen_node = weight.GetDAG().GetRoot();
    auto sample = AddMATConversion(weight.SampleTree({}, chosen_node));
    MAT::Tree mat;
    sample.View().BuildMAT(mat);
    sample.View().GetRoot().Validate(true);
    check_edge_mutations(sample.View().Const());
    Test_Move_Found_Callback callback{sample.View(), merge};
    optimized_dags.push_back(
        optimize_dag_direct(sample.View(), callback, callback, callback));
    optimized_dags.back().first.View().RecomputeCompactGenomes();
    merge.AddDAG(optimized_dags.back().first.View(), chosen_node);
  }
}

template <typename SampleDAG>
struct Single_Move_Callback_With_Hypothetical_Tree : public Move_Found_Callback {
  Single_Move_Callback_With_Hypothetical_Tree(Merge<MADAG>& merge, SampleDAG sample)
      : merge_{merge}, sample_{sample}, approved_a_move_{false} {}

  bool operator()(Profitable_Moves& move, int /*best_score_change*/,
                  std::vector<Node_With_Major_Allele_Set_Change>&
                      nodes_with_major_allele_set_change) override {
    if (!approved_a_move_) {
      // apply move to merge object.

      auto storage = [this](std::string ref_seq) {
        std::unique_lock<std::mutex> lock{mutex_};
        Assert(sample_mat_ != nullptr);
        using Storage = ExtendDAGStorage<
            DefaultDAGStorage, Extend::Nodes<Deduplicate<CompactGenome>, SampleId>,
            Extend::Edges<EdgeMutations>, Extend::DAG<ReferenceSequence>>;
        auto mat_conv = AddMATConversion(Storage{});
        mat_conv.View().BuildFromMAT(*sample_mat_, ref_seq);
        check_edge_mutations(mat_conv.View());
        mat_conv.View().RecomputeCompactGenomes();
        return SPRStorage(std::move(mat_conv));
      }(sample_.GetReferenceSequence());
      auto spr = storage.View();

      // ** create hypothetical tree
      spr.InitHypotheticalTree(move, nodes_with_major_allele_set_change);

      // ** build fragment
      auto spr_fragment = spr.GetFragment();

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

  void OnReassignedStates(const MAT::Tree&) {}

  Merge<MADAG>& merge_;
  SampleDAG sample_;
  std::mutex mutex_;
  MAT::Tree* sample_mat_ = nullptr;
  bool approved_a_move_;
};

static void test_optimizing_with_hypothetical_tree(
    const MADAGStorage& tree_shaped_dag) {
  // tbb::global_control c(tbb::global_control::max_allowed_parallelism, 1);
  // this test takes a tree and uses matOptimize to apply a single move.

  Merge<MADAG> dag_altered_in_callback{tree_shaped_dag.View().GetReferenceSequence()};
  dag_altered_in_callback.AddDAG(tree_shaped_dag.View());
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

  optimized_dag.View().RecomputeCompactGenomes();
  Merge<MADAG> two_tree_dag{tree_shaped_dag.View().GetReferenceSequence()};
  two_tree_dag.AddDAG(sample.View());
  two_tree_dag.AddDAG(optimized_dag.View());

  // check topologies of the two DAGs match in feature count
  Assert(two_tree_dag.GetResult().GetNodesCount() ==
         dag_altered_in_callback.GetResult().GetNodesCount());
  Assert(two_tree_dag.GetResult().GetEdgesCount() ==
         dag_altered_in_callback.GetResult().GetEdgesCount());
}

static auto MakeSampleDAG() {
  MADAGStorage input_storage;
  auto dag = input_storage.View();

  dag.SetReferenceSequence("GAA");

  dag.InitializeNodes(11);

  dag.AddEdge({0}, {0}, {10}, {0});
  dag.AddEdge({1}, {7}, {1}, {0}).GetMutableEdgeMutations()[{1}] = {'T', 'A'};
  dag.AddEdge({2}, {7}, {2}, {1}).GetMutableEdgeMutations()[{1}] = {'T', 'G'};
  dag.AddEdge({3}, {8}, {3}, {0}).GetMutableEdgeMutations()[{1}] = {'C', 'A'};
  dag.AddEdge({4}, {8}, {4}, {1}).GetMutableEdgeMutations()[{1}] = {'C', 'A'};
  dag.AddEdge({5}, {9}, {5}, {0}).GetMutableEdgeMutations()[{1}] = {'A', 'C'};
  dag.AddEdge({6}, {9}, {6}, {1}).GetMutableEdgeMutations()[{1}] = {'A', 'T'};
  dag.AddEdge({7}, {8}, {7}, {2}).GetMutableEdgeMutations()[{1}] = {'C', 'T'};
  dag.AddEdge({8}, {10}, {8}, {0}).GetMutableEdgeMutations()[{1}] = {'G', 'C'};
  dag.AddEdge({9}, {10}, {9}, {1}).GetMutableEdgeMutations()[{1}] = {'G', 'A'};

  dag.BuildConnections();

  dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{2}] = {'G', 'C'};
  dag.Get(EdgeId{2}).GetMutableEdgeMutations()[{2}] = {'G', 'T'};
  dag.Get(EdgeId{3}).GetMutableEdgeMutations()[{2}] = {'T', 'G'};
  dag.Get(EdgeId{4}).GetMutableEdgeMutations()[{2}] = {'T', 'C'};
  dag.Get(EdgeId{5}).GetMutableEdgeMutations()[{2}] = {'G', 'T'};
  dag.Get(EdgeId{6}).GetMutableEdgeMutations()[{2}] = {'G', 'C'};
  dag.Get(EdgeId{7}).GetMutableEdgeMutations()[{2}] = {'T', 'G'};
  dag.Get(EdgeId{8}).GetMutableEdgeMutations()[{2}] = {'A', 'T'};
  dag.Get(EdgeId{9}).GetMutableEdgeMutations()[{2}] = {'A', 'G'};

  dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{3}] = {'G', 'C'};
  dag.Get(EdgeId{2}).GetMutableEdgeMutations()[{3}] = {'G', 'T'};
  dag.Get(EdgeId{3}).GetMutableEdgeMutations()[{3}] = {'T', 'G'};
  dag.Get(EdgeId{4}).GetMutableEdgeMutations()[{3}] = {'T', 'G'};
  dag.Get(EdgeId{5}).GetMutableEdgeMutations()[{3}] = {'G', 'T'};
  dag.Get(EdgeId{6}).GetMutableEdgeMutations()[{3}] = {'G', 'C'};
  dag.Get(EdgeId{7}).GetMutableEdgeMutations()[{3}] = {'T', 'G'};
  dag.Get(EdgeId{8}).GetMutableEdgeMutations()[{3}] = {'A', 'C'};
  dag.Get(EdgeId{9}).GetMutableEdgeMutations()[{3}] = {'A', 'T'};

  dag.RecomputeCompactGenomes();
  return input_storage;
}

[[maybe_unused]] static void test_sample() {
  auto input_storage = MakeSampleDAG();
  auto dag = input_storage.View();
  auto spr_storage = SPRStorage(dag);
  auto spr = spr_storage.View();

  spr.GetRoot().Validate(true);
  spr.ApplyMove({1}, {10});
  spr.GetRoot().Validate(true);

  for (auto node : spr.GetNodes()) {
    if (not node.IsOverlaid<CompactGenome>()) {
      node.SetOverlay<CompactGenome>();
    }
  }

  spr.RecomputeCompactGenomes();
}

[[maybe_unused]] static const auto test_added0 = add_test(
    {[] {
       test_spr(Load("data/seedtree/seedtree.pb.gz", "data/seedtree/refseq.txt.gz"), 3);
     },
     "SPR: seedtree"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] {
                test_spr(Load("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
                              "data/20D_from_fasta/refseq.txt.gz"),
                         3);
              },
              "SPR: tree 20D_from_fasta"});

[[maybe_unused]] static const auto test_added2 =
    add_test({[] { test_spr(MakeSampleDAG(), 3000); }, "SPR: sample"});

[[maybe_unused]] static const auto test_added3 =
    add_test({[] { test_sample(); }, "SPR: move"});

[[maybe_unused]] static const auto test_added4 =
    add_test({[] {
                test_optimizing_with_hypothetical_tree(
                    Load("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
                         "data/20D_from_fasta/refseq.txt.gz"));
              },
              "SPR: single"});

[[maybe_unused]] static const auto test_added5 =
    add_test({[] {
                test_optimizing_with_hypothetical_tree(
                    Load("data/seedtree/seedtree.pb.gz", "data/seedtree/refseq.txt.gz")

                );
              },
              "SPR: seedtree single move"});

[[maybe_unused]] static const auto test_added6 =
    add_test({test_sample, "SPR: sample move"});

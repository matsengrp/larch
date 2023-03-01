#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/merge/merge.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/spr/spr_view.hpp"

#include <tbb/global_control.h>

template <typename DAG>
struct Test_Move_Found_Callback : public Move_Found_Callback {
  Test_Move_Found_Callback(DAG sample_dag) : sample_dag_{sample_dag} {}

  bool operator()(Profitable_Moves& move, int best_score_change,
                  std::vector<Node_With_Major_Allele_Set_Change>&
                      nodes_with_major_allele_set_change) override {
    Assert(move.src != nullptr);
    Assert(move.dst != nullptr);

    auto spr_storage = SPRStorage(sample_dag_);
    auto spr = spr_storage.View();

    Assert(sample_mat_ != nullptr);
    spr.InitHypotheticalTree(move, nodes_with_major_allele_set_change);

    for (auto node : spr.Const().GetNodes()) {
      if (node.IsRoot() or node.IsMoveNew()) {
        continue;
      }
      // std::ignore = node.Const().ComputeNewCompactGenome();
      // std::ignore = spr.GetFragment();
    }

    return move.score_change < best_score_change;
  }

  void operator()(const MAT::Tree& tree) { sample_mat_ = std::addressof(tree); }

  DAG sample_dag_;
  const MAT::Tree* sample_mat_ = nullptr;
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
  // tbb::global_control c(tbb::global_control::max_allowed_parallelism, 1);
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
    std::cout << "Sample nodes count: " << sample.GetNodesCount() << "\n";
    check_edge_mutations(sample.View());
    Test_Move_Found_Callback callback{sample.View()};
    optimized_dags.push_back(optimize_dag_direct(sample.View(), callback, callback));
    optimized_dags.back().first.View().RecomputeCompactGenomes();
    merge.AddDAG(optimized_dags.back().first.View(), chosen_node);
  }
}

static auto MakeSampleDAG() {
  MADAGStorage input_storage;
  auto dag = input_storage.View();

  dag.SetReferenceSequence("GAA");

  dag.InitializeNodes(11);

  dag.AddEdge({0}, {0}, {10}, {0});
  dag.AddEdge({1}, {7}, {1}, {0}).GetMutableEdgeMutations()[{1}] = {'T', 'A'};
  dag.AddEdge({2}, {7}, {2}, {0}).GetMutableEdgeMutations()[{1}] = {'T', 'G'};
  dag.AddEdge({3}, {8}, {3}, {0}).GetMutableEdgeMutations()[{1}] = {'C', 'A'};
  dag.AddEdge({4}, {8}, {4}, {0}).GetMutableEdgeMutations()[{1}] = {'C', 'A'};
  dag.AddEdge({5}, {9}, {5}, {0}).GetMutableEdgeMutations()[{1}] = {'A', 'C'};
  dag.AddEdge({6}, {9}, {6}, {0}).GetMutableEdgeMutations()[{1}] = {'A', 'T'};
  dag.AddEdge({7}, {8}, {7}, {0}).GetMutableEdgeMutations()[{1}] = {'C', 'T'};
  dag.AddEdge({8}, {10}, {8}, {0}).GetMutableEdgeMutations()[{1}] = {'G', 'C'};
  dag.AddEdge({9}, {10}, {9}, {0}).GetMutableEdgeMutations()[{1}] = {'G', 'A'};

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

  spr.ApplyMove({7}, {9});

  for (auto node : spr.GetNodes()) {
    if (not node.IsOverlaid<CompactGenome>()) {
      node.SetOverlay<CompactGenome>();
    }
  }

  spr.RecomputeCompactGenomes();
  
  MADAGToDOT(spr, std::cout);
}

[[maybe_unused]] static const auto test_added1 =
    add_test({[] {
                test_spr(Load("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
                              "data/20D_from_fasta/refseq.txt.gz"),
                         3);
              },
              "SPR: tree 20D_from_fasta"});

[[maybe_unused]] static const auto test_added2 =
    add_test({[] { test_spr(MakeSampleDAG(), 3); }, "SPR: sample"});

[[maybe_unused]] static const auto test_added3 =
    add_test({[] { test_sample(); }, "SPR: move"});

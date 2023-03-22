#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/merge/merge.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/spr/spr_view.hpp"

#include <fstream>
#include <tbb/global_control.h>

std::ostream& operator<<(std::ostream& os, EdgeId edge_id) {
  os << "EdgeId::" << edge_id.value;
  return os;
}

std::ostream& operator<<(std::ostream& os, NodeId node_id) {
  os << "NodeId::" << node_id.value;
  return os;
}

template <typename DAG>
struct Test_Move_Found_Callback : public Move_Found_Callback {
  Test_Move_Found_Callback(DAG sample_dag) : sample_dag_{sample_dag} {}

  bool operator()(Profitable_Moves& move, int best_score_change,
                  [[maybe_unused]] std::vector<Node_With_Major_Allele_Set_Change>&
                      nodes_with_major_allele_set_change) override {
    Assert(move.src != nullptr);
    Assert(move.dst != nullptr);
    auto storage = [this](std::string ref_seq) {
      MAT::Tree* mat = sample_mat_.load();
      Assert(mat != nullptr);
      using Storage = ExtendDAGStorage<
          DefaultDAGStorage, Extend::Nodes<Deduplicate<CompactGenome>, SampleId>,
          Extend::Edges<EdgeMutations>, Extend::DAG<ReferenceSequence>>;
      auto mat_conv = AddMATConversion(Storage{});
      mat_conv.View().BuildFromMAT(*mat, ref_seq);
      check_edge_mutations(mat_conv.View());
      mat_conv.View().RecomputeCompactGenomes();
      return SPRStorage(std::move(mat_conv));
    }(sample_dag_.GetReferenceSequence());

    auto spr = storage.View();
    spr.GetRoot().Validate(true);
    spr.InitHypotheticalTree(move, nodes_with_major_allele_set_change);
    spr.GetRoot().Validate(true);
    std::ignore = spr.GetFragment();
    return move.score_change < best_score_change;
  }

  void operator()(MAT::Tree& tree) { sample_mat_.store(std::addressof(tree)); }

  DAG sample_dag_;
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
    sample.View().GetRoot().Validate(true);
    check_edge_mutations(sample.View());
    Test_Move_Found_Callback callback{sample.View()};
    optimized_dags.push_back(optimize_dag_direct(sample.View(), callback, callback));
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
  auto [optimized_dag, optimized_mat] =
      optimize_dag_direct(sample.View(), single_move_callback, single_move_callback);

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

  MADAGToDOT(spr, std::cout);

  spr.ApplyMove({1}, {10});

  for (auto node : spr.GetNodes()) {
    if (not node.IsOverlaid<CompactGenome>()) {
      node.SetOverlay<CompactGenome>();
    }
  }

  spr.RecomputeCompactGenomes();

  MADAGToDOT(spr, std::cout);
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
                test_optimizing_with_hypothetical_tree(Load(
                    "data/seedtree/seedtree.pb.gz", "data/seedtree/refseq.txt.gz"));
              },
              "SPR: seedtree single move"});

[[maybe_unused]] static const auto test_added6 =
    add_test({test_sample, "SPR: sample move"});

// ** NEW TEST

template <typename DAG, typename SPR>
[[maybe_unused]] static void compare_dag_vs_spr_nodes(DAG dag, SPR spr,
                                                      NodeId src_node_id,
                                                      NodeId dest_node_id) {
  auto src_node = dag.Get(src_node_id);
  auto dest_node = dag.Get(dest_node_id);
  auto new_node_id = NodeId{dag.GetNodesCount()};
  for (auto node : spr.GetNodes()) {
    std::cout << "==> node: " << node.GetId() << std::endl;
    std::cout << "==> node_parents: " << node.ParentsToString() << std::endl;
    std::cout << "==> node_children: " << node.ChildrenToString() << std::endl;
    if (node.GetId().value != dag.GetNodesCount()) {
      auto dag_node = dag.Get(node.GetId());
      std::cout << "==> dag_node_parents: " << dag_node.ParentsToString() << std::endl;
      std::cout << "==> dag_node_children: " << dag_node.ChildrenToString()
                << std::endl;
    }

    // Special case: new node.
    if (node.GetId() == new_node_id) {
      std::cout << "new_node: " << node.GetId().value << std::endl;
      // Check its parens are the same as the previous parent of the dest_node.
      for (auto parent_node : dest_node.GetParentNodes()) {
        assert_true(
            node.ContainsParent(parent_node.GetId()),
            "new_node's parent is not the dest_node's previous parent after SPR.");
      }
      // Check its children are src_node and dest_node.
      for (auto child_node : {src_node, dest_node}) {
        assert_true(node.ContainsChild(child_node.GetId()),
                    "new_node's children are not src_node and dest_node after SPR.");
      }
      continue;
    }

    // Node from DAG.
    auto dag_node = dag.Get(node.GetId());
    // Cases: node can match multiple types.
    bool node_is_src = (node.GetId() == src_node.GetId());
    bool node_is_src_parent = src_node.ContainsParent(node.GetId());
    bool node_is_dest = node.GetId() == dest_node.GetId();
    bool node_is_dest_parent = dest_node.ContainsParent(node.GetId());
    bool node_is_old =
        !(node_is_src || node_is_dest || node_is_src_parent || node_is_dest_parent);

    if (node_is_old) {
      std::cout << "old_node: " << node.GetId().value << std::endl;
      // Check its parents are the same.
      for (auto parent : node.GetParentNodes()) {
        assert_true(dag_node.ContainsParent(parent.GetId()),
                    "old_node's parents are not the same after SPR.");
      }
      // Check its children are the same.
      for (auto child : node.GetChildNodes()) {
        assert_true(dag_node.ContainsChild(child.GetId()),
                    "old_node's children are not the same after SPR.");
      }
    }
    if (node_is_src) {
      std::cout << "src_node: " << node.GetId().value << std::endl;
      // Check its parent is a new node.
      for (auto parent_node : node.GetParentNodes()) {
        assert_true((parent_node.GetId() == new_node_id),
                    "src_node's parent is not new node after SPR.");
      }
      // Check its children are the same.
      for (auto child_node : node.GetChildNodes()) {
        assert_true(dag_node.ContainsChild(child_node.GetId()),
                    "src_node's children are not the same after SPR.");
      }
    } else if (node_is_dest) {
      std::cout << "dest_node: " << node.GetId().value << std::endl;
      // Check its parent is a new node.
      for (auto parent_node : node.GetParentNodes()) {
        assert_true((parent_node.GetId() == new_node_id),
                    "dest_node's parent is not new node after SPR.");
      }
      // Check its children are the same.
      for (auto child_node : node.GetChildNodes()) {
        assert_true(dag_node.ContainsChild(child_node.GetId()),
                    "dest_node's children are not the same after SPR.");
      }
    }
    if (node_is_src_parent && node_is_dest_parent) {
      std::cout << "src_and_dest_parent: " << node.GetId().value << std::endl;
      // Check its parents are the same.
      for (auto parent_node : node.GetParentNodes()) {
        test_true(dag_node.ContainsParent(parent_node.GetId()),
                  "src_parent's parents are not the same after SPR.");
      }
      // Check it has all previous children, minus src_node and dest_node, plus new
      // node.
      for (auto child_node : node.GetChildNodes()) {
        test_true(dag_node.ContainsChild(child_node.GetId()),
                  "dest_parent's children are not the same after SPR.");
      }
      test_false(node.ContainsChild(src_node.GetId()),
                 "src_parent's children incorrectly contains src_node after SPR.");
      test_false(node.ContainsChild(dest_node.GetId()),
                 "dest_parent's children incorrectly contains dest_node after SPR.");
      test_true(node.ContainsChild(new_node_id),
                "dest_parent's missing new_node after SPR.");
    } else if (node_is_src_parent) {
      std::cout << "src_parent: " << node.GetId().value << std::endl;
      if (!node_is_dest) {
        // Check its parents are the same.
        for (auto parent_node : node.GetParentNodes()) {
          test_true(dag_node.ContainsParent(parent_node.GetId()),
                    "src_parent's parents are not the same after SPR.");
        }
      }
      // Check it has all previous children, minus src_node.
      for (auto child_node : node.GetChildNodes()) {
        test_true(dag_node.ContainsChild(child_node.GetId()),
                  "src_parent's children are not same after SPR.");
      }
      test_false(node.ContainsChild(src_node.GetId()),
                 "src_parent's children incorrectly contains src_node after SPR.");
    } else if (node_is_dest_parent) {
      std::cout << "dest_parent: " << node.GetId().value << std::endl;
      if (!node_is_src) {
        // Check its parents are the same.
        for (auto parent_node : node.GetParentNodes()) {
          test_true(dag_node.ContainsParent(parent_node.GetId()),
                    "dest_parent's parents are not the same after SPR.");
        }
      }
      // Check it has all previous children, minus dest_node, plus new node.
      for (auto child_node : node.GetChildNodes()) {
        test_true(dag_node.ContainsChild(child_node.GetId()),
                  "dest_parent's children are not the same after SPR.");
      }
      test_false(node.ContainsChild(dest_node.GetId()),
                 "dest_parent's children incorrectly contains dest_node after SPR.");
      test_true(node.ContainsChild(new_node_id),
                "dest_parent's missing new_node after SPR.");
    }
  }
}

// Test that all elligible SPR moves result in valid tree states.
[[maybe_unused]] static void validate_dag_after_spr() {
  std::cout << std::endl;
  std::cout << "==> VALIDATE DAG AFTER SPR [begin] <==" << std::endl;

  auto input_storage = MakeSampleDAG();
  auto dag = input_storage.View();
  MADAGToDOT(dag, std::cout);

  size_t spr_count = 0;
  std::ofstream os;
  std::string output_folder = "_ignore/";
  std::string output_prefix = "validate_spr.";
  std::string output_ext = ".dot";
  std::string output_filename;
  // Test all valid SPR moves.
  for (auto src_node : dag.GetNodes()) {
    // Check that src_node is elligible for SPR move.
    if (src_node.IsRoot()) {
      continue;
    }
    if ((src_node.GetParentsCount() == 1) && src_node.GetSingleParent().IsRoot()) {
      continue;
    }
    for (auto dest_node : dag.GetNodes()) {
      // Check that dest_node is elligible for SPR move.
      if (dest_node.IsRoot() || (src_node.GetId() == dest_node.GetId())) {
        continue;
      }
      std::cout << "SRC: NodeId::" << src_node.GetId().value
                << ", DEST: NodeId::" << dest_node.GetId().value << std::endl;

      // Apply SPR
      auto spr_storage = SPRStorage(dag);
      auto spr = spr_storage.View();
      spr.ApplyMove(src_node.GetId(), dest_node.GetId());

      output_filename =
          output_folder + output_prefix + std::to_string(src_node.GetId().value) + "_" +
          std::to_string(dest_node.GetId().value) + ".pre_compact" + output_ext;
      std::cout << ">> WRITE DOTFILE [pre]: " << output_filename << std::endl;
      os.open(output_filename);
      MADAGToDOT(spr, os);
      os.close();

      // Update Compact Genomes.
      std::cout << ">> OVERLAY_NODES" << std::endl;
      for (auto node : spr.GetNodes()) {
        if (not node.IsOverlaid<CompactGenome>()) {
          node.SetOverlay<CompactGenome>();
        }
      }
      std::cout << ">> RECOMPUTE_NODES" << std::endl;
      // spr.RecomputeCompactGenomes();

      output_filename =
          output_folder + output_prefix + std::to_string(src_node.GetId().value) + "_" +
          std::to_string(dest_node.GetId().value) + ".post_compact" + output_ext;
      std::cout << ">> WRITE DOTFILE [post]: " << output_filename << std::endl;
      os.open(output_filename);
      MADAGToDOT(spr, os);
      os.close();

      // Test parent-child identities.
      std::cout << ">> TEST" << std::endl;
      compare_dag_vs_spr_nodes(dag, spr, src_node.GetId(), dest_node.GetId());
    }
    break;
  }

  std::cout << "==> VALIDATE DAG AFTER SPR [end] <==" << std::endl;
}

[[maybe_unused]] static const auto test_added7 =
    add_test({validate_dag_after_spr, "SPR: validate DAG after SPR"});

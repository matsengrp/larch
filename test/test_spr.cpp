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

template <typename DAG, typename SPR, typename Node>
[[maybe_unused]] static bool test_compare_dag_vs_spr_nodes(
    DAG dag, SPR spr, Node src_dag, Node dest_dag, std::vector<size_t> child_counts) {
  bool passes_test = true;
  auto src_spr = spr.Get(src_dag.GetId());
  auto dest_spr = spr.Get(dest_dag.GetId());

  // Nodes that have collapsed during move (these indices can be reused).
  std::set<NodeId> collapsed_nodes;
  // Include new node.
  auto new_node_id = NodeId{dag.GetNodesCount()};
  collapsed_nodes.insert(new_node_id);
  // Check if an edge should collapse on move. Node should collapse if it has a zero or
  // one child and has a single parent after move.
  for (auto parent_node : src_dag.GetParentNodes()) {
    if ((child_counts[parent_node.GetId().value] <= 2) &&
        (parent_node.GetParentsCount() == 1)) {
      collapsed_nodes.insert(parent_node.GetId());
    }
  }

  // Nodes that are expected to change from move.
  std::set<NodeId> changed_nodes;
  for (auto node : {src_dag, dest_dag}) {
    changed_nodes.insert(node.GetId());
    for (auto parent_node : node.GetParentNodes()) {
      changed_nodes.insert(parent_node.GetId());
    }
  }

  // Check that src_node and dest_node have the same parents.
  for (auto src_parent_node : src_spr.GetParentNodes()) {
    if (!dest_spr.ContainsParent(src_parent_node.GetId())) {
      std::cout
          << "TEST_FAILED: src_node and dest_node do not share src_parent_node -- "
          << src_parent_node.GetId().value << std::endl;
      passes_test = false;
    }
  }
  for (auto dest_parent_node : dest_spr.GetParentNodes()) {
    if (!src_spr.ContainsParent(dest_parent_node.GetId())) {
      std::cout
          << "TEST_FAILED: src_node and dest_node do not share dest_parent_node -- "
          << dest_parent_node.GetId().value << std::endl;
      passes_test = false;
    }
  }

  // Check that parents of src_node and dest_node are either previous dest_parent nodes
  // or collapsed nodes.
  for (auto parent_node : src_spr.GetParentNodes()) {
    if (!(dest_dag.ContainsParent(parent_node.GetId())) &&
        !(collapsed_nodes.find(parent_node.GetId()) != collapsed_nodes.end())) {
      std::cout << "TEST_FAILED: dest_and_src_parent is not from previous dest_node or "
                   "collapsed nodes -- "
                << parent_node.GetId().value << std::endl;
      passes_test = false;
    }
  }

  // Check that all connections are reflective (every node's child is also that node's
  // parent).
  for (auto node : spr.GetNodes()) {
    for (auto child_node : node.GetChildNodes()) {
      if (!child_node.ContainsParent(node.GetId())) {
        std::cout << "TEST_FAILED: node and child_node are not reflective -- "
                  << node.GetId().value << " " << child_node.GetId().value << std::endl;
        passes_test = false;
      }
    }
    for (auto parent_node : node.GetParentNodes()) {
      if (!parent_node.ContainsChild(node.GetId())) {
        std::cout << "TEST_FAILED: node and parent_node are not reflective -- "
                  << node.GetId().value << " " << parent_node.GetId().value
                  << std::endl;
        passes_test = false;
      }
    }
  }

  // Check that all other nodes are not modified from pre-move DAG.
  for (auto dag_node : dag.GetNodes()) {
    if ((changed_nodes.find(dag_node.GetId()) != changed_nodes.end()) ||
        (collapsed_nodes.find(dag_node.GetId()) != collapsed_nodes.end())) {
      continue;
    }
    auto spr_node = spr.Get(dag_node.GetId());
    for (auto dag_parent : dag_node.GetParentNodes()) {
      if (!spr_node.ContainsParent(dag_parent.GetId()) &&
          !(collapsed_nodes.find(dag_parent) != collapsed_nodes.end())) {
        std::cout << "TEST_FAILED: unchanged spr_node missing child from dag_node -- "
                  << dag_node.GetId().value << "->" << dag_parent.GetId().value
                  << std::endl;
        passes_test = false;
      }
    }
    for (auto dag_child : dag_node.GetChildNodes()) {
      if (!spr_node.ContainsChild(dag_child.GetId()) &&
          !(collapsed_nodes.find(dag_child) != collapsed_nodes.end())) {
        std::cout << "TEST_FAILED: unchanged spr_node missing parent from dag_node -- "
                  << dag_node.GetId().value << "->" << dag_child.GetId().value
                  << std::endl;
        passes_test = false;
      }
    }
  }

  return passes_test;
}

template <typename DAG>
std::vector<size_t> get_child_counts(DAG dag) {
  std::vector<size_t> child_counts(dag.GetNodesCount());
  for (auto node : dag.GetNodes()) {
    size_t child_count = 0;
    for (auto child : node.GetChildren()) {
      std::ignore = child;
      child_count++;
    }
    child_counts[node.GetId().value] = child_count;
  }
  return child_counts;
}

template <typename Node>
std::set<Node> get_all_parents_of_node(Node child_node) {
  std::set<Node> parent_nodes, new_nodes;
  for (auto parent_node : child_node.GetParentNodes()) {
    parent_nodes.insert(parent_node);
  }
  new_nodes = parent_nodes;
  while (!new_nodes.empty()) {
    std::set<Node> new_new_nodes;
    for (auto new_node : new_nodes) {
      for (auto parent_node : new_node.GetParentNodes()) {
        parent_nodes.insert(parent_node);
        new_new_nodes.insert(parent_node);
      }
    }
    new_nodes = new_new_nodes;
  }
  return parent_nodes;
}

template <typename DAG, typename Node>
bool is_valid_spr_move(DAG dag, Node src_node, Node dest_node) {
  if (src_node.GetId() == dest_node.GetId()) {
    return false;
  }
  if (src_node.IsRoot() || dest_node.IsRoot()) {
    return false;
  }
  for (auto parent_node : src_node.GetParentNodes()) {
    if (parent_node.GetId() == dest_node.GetId()) {
      return false;
    }
  }
  const auto dest_parents = get_all_parents_of_node(dest_node);
  if (dest_parents.find(src_node) != dest_parents.end()) {
    return false;
  }
  return true;
}

// Test that all elligible SPR moves result in valid tree states.
[[maybe_unused]] static void validate_dag_after_spr() {
  std::ofstream os;
  std::string output_filename;
  std::string output_folder = "_ignore/";
  std::string output_prefix = "validate_spr.";
  std::string output_ext = ".dot";

  std::cout << std::endl;
  std::cout << "==> VALIDATE DAG AFTER SPR [begin] <==" << std::endl;

  auto dag_storage = MakeSampleDAG();
  auto dag = dag_storage.View();
  auto child_counts = get_child_counts(dag);

  output_filename = output_folder + "sample_dag" + output_ext;
  std::cout << ">> WRITE DOTFILE [pre]: " << output_filename << std::endl;
  os.open(output_filename);
  MADAGToDOT(dag, os);
  os.close();

  // Test all valid SPR moves.
  for (auto src_node : dag.GetNodes()) {
    for (auto dest_node : dag.GetNodes()) {
      if (!is_valid_spr_move(dag, src_node, dest_node)) {
        std::cout << "SKIPPING..." << std::endl;
        continue;
      }

      std::cout << "SRC: NodeId::" << src_node.GetId().value
                << ", DEST: NodeId::" << dest_node.GetId().value << std::endl;

      // Apply SPR
      auto spr_storage = SPRStorage(dag);
      auto spr = spr_storage.View();
      spr.ApplyMove(src_node.GetId(), dest_node.GetId());

      // Update Compact Genomes.
      for (auto node : spr.GetNodes()) {
        if (not node.IsOverlaid<CompactGenome>()) {
          node.SetOverlay<CompactGenome>();
        }
      }
      spr.RecomputeCompactGenomes();

      output_filename = output_folder + "sample_dag." +
                        std::to_string(src_node.GetId().value) + "_" +
                        std::to_string(dest_node.GetId().value) + output_ext;
      std::cout << ">> WRITE DOTFILE [post]: " << output_filename << std::endl;
      os.open(output_filename);
      MADAGToDOT(spr, os);
      os.close();

      // Test parent-child identities.
      test_compare_dag_vs_spr_nodes(dag, spr, src_node, dest_node, child_counts);
    }
  }
  std::cout << "==> VALIDATE DAG AFTER SPR [end] <==" << std::endl;
}

[[maybe_unused]] static const auto test_added7 =
    add_test({validate_dag_after_spr, "SPR: validate DAG after SPR"});

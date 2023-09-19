#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/merge/merge.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/spr/spr_view.hpp"
#include "larch/spr/lca.hpp"
#include "sample_dag.hpp"

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
      mat_conv.View().RecomputeCompactGenomes(true);
      return SPRStorage(std::move(mat_conv));
    }(sample_dag_.GetReferenceSequence());

    auto spr = storage.View();
    spr.GetRoot().Validate(true);

    if (spr.InitHypotheticalTree(move, nodes_with_major_allele_set_change)) {
      spr.GetRoot().Validate(true);

      auto fragment = spr.GetFragment();

      std::scoped_lock<std::mutex> lock{merge_mtx_};
      merge_.AddFragment(spr, fragment.first, fragment.second);
    } else {
      return false;
    }
    return move.score_change < best_score_change;
  }

  void operator()(MAT::Tree& tree) {
    decltype(AddMATConversion(Storage{{}})) storage;
    storage.View().BuildFromMAT(tree, sample_dag_.GetReferenceSequence());
    storage.View().RecomputeCompactGenomes(true);
    {
      std::scoped_lock<std::mutex> lock{merge_mtx_};
      merge_.AddDAG(storage.View());
      sample_mat_.store(std::addressof(tree));
      merge_.ComputeResultEdgeMutations();
    }
    // StoreDAGToProtobuf(merge_.GetResult(), "radius_iter.pb");
  }

  void OnReassignedStates(MAT::Tree& tree) {
    reassigned_states_storage_.View().BuildFromMAT(tree,
                                                   sample_dag_.GetReferenceSequence());
    check_edge_mutations(reassigned_states_storage_.View().Const());
    reassigned_states_storage_.View().RecomputeCompactGenomes(true);
    {
      std::scoped_lock<std::mutex> lock{merge_mtx_};
      merge_.AddDAG(reassigned_states_storage_.View());
      merge_.ComputeResultEdgeMutations();
    }
  }

  DAG sample_dag_;
  MergeT& merge_;
  decltype(AddMATConversion(Storage{{}})) reassigned_states_storage_ =
      AddMATConversion(Storage{});
  std::atomic<MAT::Tree*> sample_mat_ = nullptr;
  std::mutex merge_mtx_;
};

[[maybe_unused]] static MADAGStorage Load(std::string_view input_dag_path,
                                          std::string_view refseq_path) {
  std::string reference_sequence = LoadReferenceSequence(refseq_path);
  MADAGStorage input_dag_storage =
      LoadTreeFromProtobuf(input_dag_path, reference_sequence);
  input_dag_storage.View().RecomputeCompactGenomes(true);
  return input_dag_storage;
}
// Use the base DAG and the src and dest nodes.
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
  // or collapsed nodes (collapsed node ids can be repurposed for new node ids).
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

  // Check that all other nodes are not modified from pre-move DAG. Allows for the
  // exception when parent or child is a collapsed node.
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

template <typename Node>
bool is_valid_spr_move(Node src_node, Node dest_node) {
  if (src_node.GetId() == dest_node.GetId()) {
    return false;
  }
  if (src_node.IsUA() || dest_node.IsUA()) {
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
  const auto src_parents = get_all_parents_of_node(src_node);
  const auto dest_parent_node = dest_node.GetSingleParent().GetParent();
  if (src_parents.find(dest_parent_node) != src_parents.end()) {
    return false;
  }
  return true;
}

// Test that all elligible SPR moves result in valid tree states.
[[maybe_unused]] static void validate_dag_after_spr(const std::string& dag_name,
                                                    bool write_dot_files = false) {
  std::ofstream os;
  std::string output_filename;
  std::string output_folder = "_ignore/";
  std::string output_ext = ".dot";

  auto dag_storage = AddMATConversion(MakeSampleDAG());
  auto dag = dag_storage.View();
  MAT::Tree tree;
  dag.BuildMAT(tree);
  auto child_counts = get_child_counts(dag);

  if (write_dot_files) {
    output_filename = output_folder + dag_name + output_ext;
    std::cout << ">> WRITE DOTFILE: " << output_filename << std::endl;
    os.open(output_filename);
    MADAGToDOT(dag, os);
    os.close();
  }

  // Test all valid SPR moves.
  for (auto src_node : dag.GetNodes()) {
    for (auto dest_node : dag.GetNodes()) {
      if (!is_valid_spr_move(src_node, dest_node)) {
        continue;
      }

      // Apply SPR
      auto spr_storage = SPRStorage(dag);
      auto spr = spr_storage.View();
      spr.GetRoot().Validate(true);
      LCA lca = FindLCA(src_node, dest_node);
      auto move_result = spr.ApplyMove(lca.lca, src_node.GetId(), dest_node.GetId());
      if (move_result.first.value != NoId) {
        // Update Compact Genomes.
        for (auto node : spr.GetNodes()) {
          if (not node.IsOverlaid<CompactGenome>()) {
            node.SetOverlay<CompactGenome>();
          }
        }
        spr.RecomputeCompactGenomes(true);

        if (write_dot_files) {
          std::cout << "SRC: NodeId::" << src_node.GetId().value
                    << ", DEST: NodeId::" << dest_node.GetId().value << std::endl;
          output_filename = output_folder + dag_name + "." +
                            std::to_string(src_node.GetId().value) + "_" +
                            std::to_string(dest_node.GetId().value) + output_ext;
          std::cout << ">> WRITE DOTFILE [post]: " << output_filename << std::endl;
          os.open(output_filename);
          MADAGToDOT(spr, os);
          os.close();
        }

        assert_true(
            test_compare_dag_vs_spr_nodes(dag, spr, src_node, dest_node, child_counts),
            "DAG '" + dag_name + "' created invalid DAG after move " +
                std::to_string(src_node.GetId().value) + " -> " +
                std::to_string(dest_node.GetId().value));
      }
    }
  }
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[]() { validate_dag_after_spr("sample_dag", true); },
              "SPR: validate DAG after SPR (sample_dag)"});

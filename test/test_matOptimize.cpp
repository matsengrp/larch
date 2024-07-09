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
  auto GetMATNodeToCGMap() { return std::map<MAT::Node*, CompactGenome>{}; }
};

[[maybe_unused]] static auto choose_root = [](const auto& subtree_weight) {
  return subtree_weight.GetDAG().GetRoot();
};

[[maybe_unused]] static auto choose_random = [](const auto& weight) {
  std::random_device random_device;
  std::mt19937 random_generator(random_device());

  NodeId node_id;
  do {
    TestAssert(weight.GetDAG().GetNodesCount() > 0);
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

std::vector<std::vector<const SampleId*>> clades_union(
    const std::vector<std::vector<const SampleId*>>& lhs,
    const std::vector<std::vector<const SampleId*>>& rhs) {
  std::vector<std::vector<const SampleId*>> result;

  for (auto [lhs_clade, rhs_clade] : ranges::views::zip(lhs, rhs)) {
    std::vector<const SampleId*> clade{lhs_clade};
    clade.insert(clade.end(), rhs_clade.begin(), rhs_clade.end());
    ranges::sort(clade);
    ranges::unique(clade);
    result.push_back(std::move(clade));
  }

  ranges::sort(result);
  return result;
}

std::vector<std::vector<const SampleId*>> clades_difference(
    const std::vector<std::vector<const SampleId*>>& lhs,
    const std::vector<std::vector<const SampleId*>>& rhs) {
  std::vector<std::vector<const SampleId*>> result;

  for (auto [lhs_clade, rhs_clade] : ranges::views::zip(lhs, rhs)) {
    std::vector<const SampleId*> clade;
    std::set_difference(lhs_clade.begin(), lhs_clade.end(), rhs_clade.begin(),
                        rhs_clade.end(), std::inserter(clade, clade.begin()));
    ranges::sort(clade);
    ranges::unique(clade);
    result.push_back(std::move(clade));
  }

  ranges::sort(result);
  return result;
}

// NOLINTNEXTLINE(cppcoreguidelines-virtual-class-destructor)
template <typename SampleDAG>
struct Larch_Move_Found_Callback : public Move_Found_Callback {
  Larch_Move_Found_Callback(const Merge& merge, SampleDAG sample)
      : merge_{merge}, sample_{sample}, move_score_coeffs_{1, 1} {}
  Larch_Move_Found_Callback(
      const Merge& merge, SampleDAG sample,
      std::pair<int, int> move_score_coeffs)  // NOLINT(modernize-pass-by-value)
      : merge_{merge}, sample_{sample}, move_score_coeffs_{move_score_coeffs} {}
  bool operator()(Profitable_Moves& move, int /* best_score_change */,
                  std::vector<Node_With_Major_Allele_Set_Change>&
                  /* node_with_major_allele_set_change */) override {
    int node_id_map_count = 0;
    if (move_score_coeffs_.first != 0) {
      NodeId src_id = ToMergedNodeId(move.src);
      NodeId dst_id = ToMergedNodeId(move.dst);
      NodeId lca_id = ToMergedNodeId(move.LCA);

      const auto& src_clades =
          merge_.GetResultNodeLabels().at(src_id).GetLeafSet()->GetClades();
      const auto& dst_clades =
          merge_.GetResultNodeLabels().at(dst_id).GetLeafSet()->GetClades();

      MAT::Node* curr_node = move.src;
      while (not(curr_node->node_id == lca_id.value)) {
        MergeDAG::NodeView node = merge_.GetResult().Get(NodeId{0 /*FIXME*/});
        const auto& clades =
            merge_.GetResultNodeLabels().at(node).GetLeafSet()->GetClades();
        if (not merge_.ContainsLeafset(clades_difference(clades, src_clades))) {
          ++node_id_map_count;
        }
        curr_node = curr_node->parent;
        if (curr_node == nullptr) {
          break;
        }
      }

      curr_node = move.dst;
      while (not(curr_node->node_id == lca_id.value)) {
        MergeDAG::NodeView node = merge_.GetResult().Get(NodeId{0 /*FIXME*/});
        const auto& clades =
            merge_.GetResultNodeLabels().at(node).GetLeafSet()->GetClades();
        if (not merge_.ContainsLeafset(clades_union(clades, dst_clades))) {
          ++node_id_map_count;
        }
        curr_node = curr_node->parent;
        if (curr_node == nullptr) {
          break;
        }
      }
    }

    move.score_change = move_score_coeffs_.second * move.score_change -
                        move_score_coeffs_.first * node_id_map_count;
    return move.score_change <= 0;
  }

  void MergeNodeIDs(std::map<MATNodePtr, NodeId>&& node_id_map) {
    node_id_map_.merge(std::forward<decltype(node_id_map)>(node_id_map));
  }

  void OnReassignedStates(const MAT::Tree& tree) {
    for (auto leaf_node : tree.get_leaves()) {
      auto new_cg = sample_.GetNodeFromMAT(leaf_node).GetCompactGenome().Copy();
      mat_node_to_cg_map_[leaf_node] = new_cg.Copy();
    }
  }

  auto GetMATNodeToCGMap() { return std::map<MAT::Node*, CompactGenome>{}; }

  std::map<MAT::Node*, CompactGenome> mat_node_to_cg_map_;

 private:
  NodeId ToMergedNodeId(MATNodePtr node) {
    auto it = node_id_map_.find(node);
    if (it != node_id_map_.end()) {
      return it->second;
    }
    return sample_.GetNodeFromMAT(node).GetOriginalId();
  }

  const Merge& merge_;
  SampleDAG sample_;
  const std::pair<int, int> move_score_coeffs_;
  std::map<MATNodePtr, NodeId> node_id_map_;
};

template <typename ChooseNode>
static void test_matOptimize(std::string_view input_dag_path,
                             std::string_view refseq_path, size_t count,
                             ChooseNode& choose_node) {
  std::string reference_sequence = LoadReferenceSequence(refseq_path);
  MADAGStorage<> input_dag_storage =
      LoadTreeFromProtobuf(input_dag_path, reference_sequence);
  input_dag_storage.View().RecomputeCompactGenomes(true);
  MADAG input_dag = input_dag_storage.View();
  Merge merge{input_dag.GetReferenceSequence()};
  merge.AddDAGs(std::vector{input_dag});
  std::vector<decltype(AddMappedNodes(
      AddMATConversion(MADAGStorage<>::EmptyDefault())))>
      optimized_dags;

  for (size_t i = 0; i < count; ++i) {
    merge.ComputeResultEdgeMutations();
    SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};

    auto chosen_node = choose_node(weight);
    bool subtrees = not chosen_node.IsUA();
    auto sample = AddMATConversion(weight.SampleTree({}, chosen_node));
    MAT::Tree mat;
    sample.View().BuildMAT(mat);
    std::cout << "Sample nodes count: " << sample.View().GetNodesCount() << "\n";
    check_edge_mutations(sample.View());
    int move_coeff_nodes = 1;
    int move_coeff_pscore = 1;
    Larch_Move_Found_Callback callback{
        merge, sample.View(), {move_coeff_nodes, move_coeff_pscore}};
    /* StoreTreeToProtobuf(sample.View(), "before_optimize_dag.pb"); */
    auto radius_callback = [&](MAT::Tree& tree) -> void {
      auto temp_result =
          AddMappedNodes(AddMATConversion(MADAGStorage<>::EmptyDefault()));
      temp_result.View().BuildFromMAT(tree, merge.GetResult().GetReferenceSequence());
      temp_result.View().RecomputeCompactGenomes(true);
      optimized_dags.push_back(std::move(temp_result));
      auto result = optimized_dags.back().View();
      std::map<MATNodePtr, NodeId> full_map = [&] {
        if (subtrees) {
          merge.AddDAGs(std::vector{result}, merge.GetResult().Get(chosen_node));
        } else {
          merge.AddDAGs(std::vector{result});
        }
        // mat_node_map is not the identity, so all pairs in mat_node_map must be used
        // to build remaped
        std::map<MATNodePtr, NodeId> remaped;
        for (auto node : result.GetNodes()) {
          if (node.IsUA()) {
            continue;
          }
          remaped.insert({node.GetMATNode(), node.GetOriginalId()});
        }
        return remaped;
      }();
      callback.MergeNodeIDs(std::move(full_map));
    };
    optimize_dag_direct(sample.View(), callback, radius_callback, callback);
  }
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] {
                test_matOptimize("data/seedtree/seedtree.pb.gz",
                                 "data/seedtree/refseq.txt.gz", 3, choose_root);
              },
              "matOptimize: seedtree"});

#pragma once

#include <mutex>
#include <shared_mutex>

#include "larch/spr/spr_view.hpp"
#include "larch/merge/merge.hpp"
#include "larch/mat_view.hpp"

/**
 * @brief Base class that adds batching functionality to Move_Found_Callback for optimized parallel processing.
 * 
 * BatchingCallback is designed to be inherited from using CRTP, derived classes must implement OnMove().
 * It enhances any Move_Found_Callback with batching capabilities, which is an optimization technique that
 * collects multiple SPR moves and merges them in parallel once a threshold is reached.
 * 
 * @tparam CRTP The derived class type for compile-time polymorphism
 */
template <typename CRTP>
class BatchingCallback : public Move_Found_Callback {
 public:
  BatchingCallback(Merge& merge);
  BatchingCallback(Merge& merge, bool collapse_empty_fragment_edges);

  virtual ~BatchingCallback() {
    auto lock = WriteLock(mat_mtx_);
#if USE_MAT_VIEW
    sample_mat_tree_.delete_nodes();
#else
    if (sample_mat_storage_ != nullptr) {
      sample_mat_storage_->View().GetMutableMAT().delete_nodes();
      sample_mat_storage_ = nullptr;
    }
#endif
  }

#if USE_MAT_VIEW
  using Storage = MergeDAGStorage<>;
  using MATStorage = UncondensedMergeDAGStorage;
  using SPRType = decltype(AddSPRStorage(std::declval<MATStorage>().View()));
  using ReassignedStatesStorage =
      decltype(AddMappedNodes(AddMATConversion(MergeDAGStorage<>::EmptyDefault())));
  using FragmentType = FragmentStorage<decltype(std::declval<SPRType>().View())>;
#else
  using Storage = MergeDAGStorage<>;
  using MATStorage = decltype(AddMATConversion(Storage::EmptyDefault()));
  using SPRType = decltype(AddSPRStorage(std::declval<MATStorage>().View()));
  using ReassignedStatesStorage =
      decltype(AddMappedNodes(AddMATConversion(Storage::EmptyDefault())));
  using FragmentType = FragmentStorage<decltype(std::declval<SPRType>().View())>;
#endif

  bool operator()(Profitable_Moves& move, int best_score_change,
                  std::vector<Node_With_Major_Allele_Set_Change>&
                      nodes_with_major_allele_set_change) override;
  void operator()(MAT::Tree& tree);

  void OnReassignedStates(MAT::Tree& tree);
  const GrowableHashMap<std::string, CompactGenome>& GetSampleIdToCGMap() const;

 protected:
  Merge& GetMerge();
  size_t GetAppliedMovesCount();
  auto GetMappedStorage();

 private:
  struct MoveStorage {
    MOVE_ONLY(MoveStorage);
    MoveStorage(std::unique_ptr<SPRType> spr_in, std::unique_ptr<FragmentType> fragment_in)
        : spr{std::move(spr_in)}, fragment{std::move(fragment_in)} {}
    std::unique_ptr<SPRType> spr;
    std::unique_ptr<FragmentType> fragment;
  };

#if USE_MAT_VIEW
  static MAT::Node* CopyNode(const MAT::Node* in, MAT::Node* parent) {
    MAT::Node* out = new MAT::Node{in->node_id};
    out->branch_length = in->branch_length;
    out->parent = parent;
    out->mutations = in->mutations;
    out->have_masked = in->have_masked;
    out->children.reserve(in->children.size());
    for (auto* c : in->children) {
      out->children.push_back(CopyNode(c, out));
    }
    return out;
  }

  static MAT::Tree CopyTree(const MAT::Tree& in) {
    MAT::Tree out;
    out.root = CopyNode(in.root, nullptr);
    out.node_names = in.node_names;
    out.node_name_to_idx_map = in.node_name_to_idx_map;
    out.node_idx = in.node_idx;
    out.num_nodes = in.num_nodes;
    out.root_ident = in.root_ident;
    out.max_level = in.max_level;
    out.condensed_nodes = in.condensed_nodes;
    out.curr_internal_node = in.curr_internal_node;
    out.all_nodes.resize(in.all_nodes.size());
    for (auto node : out.depth_first_expansion()) {
      out.all_nodes[node->node_id] = node;
    }
    return out;
  }

  void SetSample(MAT::Tree& tree, std::string ref_seq) {
    sample_mat_tree_.delete_nodes();
    sample_mat_tree_ = CopyTree(tree);
    sample_refseq_ = std::move(ref_seq);
  }
  MAT::Tree sample_mat_tree_;
  std::string sample_refseq_;

  void CreateMATViewStorage(MAT::Tree& tree, std::string_view ref_seq);
#else
  void CreateMATStorage(MAT::Tree& tree, std::string_view ref_seq);
#endif

  Merge& merge_;
  bool collapse_empty_fragment_edges_;

  std::unique_ptr<ReassignedStatesStorage> reassigned_states_storage_ =
      std::make_unique<ReassignedStatesStorage>(
          AddMappedNodes(AddMATConversion(MergeDAGStorage<>::EmptyDefault())));
  std::unique_ptr<MATStorage> sample_mat_storage_;

  std::atomic<size_t> applied_moves_count_;
  std::shared_mutex mat_mtx_;
  std::mutex merge_mtx_;
  Reduction<std::deque<MoveStorage>> moves_batch_{32};
};

#include "larch/impl/spr/batching_callback_impl.hpp"

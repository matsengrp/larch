#pragma once

#include <mutex>
#include <shared_mutex>

#include "larch/spr/spr_view.hpp"
#include "larch/merge/merge.hpp"
#include "larch/mat_view.hpp"

template <typename CRTP, typename SampleDAG>
class BatchingCallback : public Move_Found_Callback {
 public:
  BatchingCallback(Merge& merge, SampleDAG sample_dag);
  BatchingCallback(Merge& merge, SampleDAG sample_dag,
                   bool collapse_empty_fragment_edges);

  virtual ~BatchingCallback() {
    auto lock = WriteLock(mat_mtx_);
#if USE_MAT_VIEW
// TODO
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
  using SPRType = decltype(AddSPRStorage(std::declval<MATStorage>()));
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
    std::unique_ptr<SPRType> spr;
    std::unique_ptr<FragmentType> fragment;
  };

#if USE_MAT_VIEW
  void SetSample(MAT::Tree& tree, std::string ref_seq) {

  auto dag_from_mat = AddMATConversion(Storage::EmptyDefault());
  dag_from_mat.View().BuildFromMAT(tree, ref_seq);
  dag_from_mat.View().BuildMAT(sample_mat_tree_);

    //sample_mat_tree_ = tree.copy_tree();
    sample_refseq_ = std::move(ref_seq);
  }
  MAT::Tree sample_mat_tree_;
  std::string sample_refseq_;

  UncondensedMergeDAGStorage CreateMATViewStorage();
#else
  void CreateMATStorage(MAT::Tree& tree, std::string_view ref_seq);
#endif

  Merge& merge_;
  SampleDAG sample_dag_;
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

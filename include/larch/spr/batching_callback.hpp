#pragma once

#include <mutex>
#include <shared_mutex>

#include "larch/spr/spr_view.hpp"
#include "larch/merge/merge.hpp"

template <typename CRTP, typename SampleDAG>
class BatchingCallback : public Move_Found_Callback {
 public:
  BatchingCallback(Merge& merge, SampleDAG sample_dag);
  BatchingCallback(Merge& merge, SampleDAG sample_dag,
                   bool collapse_empty_fragment_edges);

  virtual ~BatchingCallback() {
    auto lock = WriteLock(mat_mtx_);
    if (sample_mat_storage_ != nullptr) {
      sample_mat_storage_->View().GetMutableMAT().delete_nodes();
      sample_mat_storage_ = nullptr;
    }
  }

  using Storage = MergeDAGStorage<>;  // TODO MADAG storage?
  using MATStorage = decltype(AddMATConversion(Storage::EmptyDefault()));
  using SPRType =
      decltype(AddSPRStorage(AddMATConversion(Storage::EmptyDefault()).View()));
  using ReassignedStatesStorage =
      decltype(AddMappedNodes(AddMATConversion(Storage::EmptyDefault())));
  using FragmentType = FragmentStorage<decltype(std::declval<SPRType>().View())>;

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

  void CreateMATStorage(MAT::Tree& tree, std::string_view ref_seq);

  Merge& merge_;
  SampleDAG sample_dag_;
  bool collapse_empty_fragment_edges_;
  std::atomic<size_t> applied_moves_count_;
  std::unique_ptr<ReassignedStatesStorage> reassigned_states_storage_ =
      std::make_unique<ReassignedStatesStorage>(
          AddMappedNodes(AddMATConversion(Storage::EmptyDefault())));
  std::shared_mutex mat_mtx_;
  std::unique_ptr<MATStorage> sample_mat_storage_;
  std::mutex merge_mtx_;
  Reduction<std::deque<MoveStorage>> moves_batch_{32};
};

#include "larch/impl/spr/batching_callback_impl.hpp"

#pragma once

#include <mutex>
#include <shared_mutex>

#include "larch/spr/spr_view.hpp"
#include "larch/merge/merge.hpp"

template <typename CRTP, typename SampleDAG>
class BatchingCallback : public Move_Found_Callback {
 public:
  BatchingCallback(Merge& merge, SampleDAG sample_dag);

  virtual ~BatchingCallback() {}

  using Storage = MergeDAGStorage;
  using MATStorage = decltype(AddMATConversion(Storage{{}}));
  using SPRType = decltype(SPRStorage(AddMATConversion(Storage{{}}).View()));

  bool operator()(Profitable_Moves& move, int best_score_change,
                  std::vector<Node_With_Major_Allele_Set_Change>&
                      nodes_with_major_allele_set_change) override;

  void operator()(MAT::Tree& tree);

  void OnReassignedStates(MAT::Tree& tree);

 protected:
  Merge& GetMerge();
  auto GetMappedStorage();

 private:
  void CreateMATStorage(MAT::Tree& tree, std::string_view ref_seq);

  Merge& merge_;
  std::decay_t<SampleDAG> sample_dag_;
  decltype(AddMappedNodes(AddMATConversion(Storage{{}}))) reassigned_states_storage_ =
      AddMappedNodes(AddMATConversion(Storage{{}}));
  std::shared_mutex mat_mtx_;
  std::unique_ptr<MATStorage> sample_mat_storage_;
  std::mutex merge_mtx_;
  tbb::concurrent_vector<SPRType> batch_storage_;
  tbb::concurrent_vector<Fragment<decltype(std::declval<SPRType>().View())>> batch_;
};

#include "larch/impl/spr/batching_callback_impl.hpp"

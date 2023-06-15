#pragma once

#include <mutex>
#include <shared_mutex>

#include "larch/spr/spr_view.hpp"
#include "larch/merge/merge.hpp"

template <typename DAG, typename MergeT>
struct BatchingCallback : public Move_Found_Callback {
  BatchingCallback(DAG sample_dag, MergeT& merge)
      : sample_dag_{sample_dag}, merge_{merge} {};

  virtual ~BatchingCallback() {}

  using Storage = MergeDAGStorage;
  using MATStorage = decltype(AddMATConversion(Storage{{}}));
  using SPRType = decltype(SPRStorage(AddMATConversion(Storage{{}}).View()));

  bool operator()(Profitable_Moves& move, int best_score_change,
                  std::vector<Node_With_Major_Allele_Set_Change>&
                      nodes_with_major_allele_set_change) override;

  void CreateMATStorage(MAT::Tree& tree, std::string_view ref_seq);

  void operator()(MAT::Tree& tree);

  void OnReassignedStates(MAT::Tree& tree);

  DAG sample_dag_;
  MergeT& merge_;
  decltype(AddMATConversion(Storage{{}})) reassigned_states_storage_ =
      AddMATConversion(Storage{{}});
  std::shared_mutex mat_mtx_;
  std::unique_ptr<MATStorage> sample_mat_storage_;
  std::mutex merge_mtx_;
  tbb::concurrent_vector<SPRType> batch_storage_;
  tbb::concurrent_vector<Fragment<decltype(std::declval<SPRType>().View())>> batch_;
};

#include "larch/impl/spr/batching_callback_impl.hpp"

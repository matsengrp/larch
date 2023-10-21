#pragma once

#include <mutex>
#include <shared_mutex>
#include <tbb/concurrent_unordered_map.h>

#include "larch/spr/spr_view.hpp"
#include "larch/merge/merge.hpp"

template <typename K, typename V>
using ConcurrentUnorderedMap =
    tbb::concurrent_unordered_map<K, V, std::hash<K>, std::equal_to<K>>;

template <typename CRTP, typename SampleDAG>
class BatchingCallback : public Move_Found_Callback {
 public:
  BatchingCallback(Merge& merge, SampleDAG sample_dag);
  BatchingCallback(Merge& merge, SampleDAG sample_dag,
                   bool collapse_empty_fragment_edges);

  virtual ~BatchingCallback() {}

  using Storage = MergeDAGStorage;
  using MATStorage = decltype(AddMATConversion(Storage{{}}));
  using SPRType = decltype(AddSPRStorage(AddMATConversion(Storage{{}}).View()));
  using ReassignedStatesStorage =
      decltype(AddMappedNodes(AddMATConversion(Storage{{}})));

  bool operator()(Profitable_Moves& move, int best_score_change,
                  std::vector<Node_With_Major_Allele_Set_Change>&
                      nodes_with_major_allele_set_change) override;

  void operator()(MAT::Tree& tree);

  void OnReassignedStates(MAT::Tree& tree);
  const ConcurrentUnorderedMap<std::string, CompactGenome>& GetSampleIdToCGMap() const;

 protected:
  Merge& GetMerge();
  ArbitraryInt GetAppliedMovesCount();
  auto GetMappedStorage();

 private:
  void CreateMATStorage(MAT::Tree& tree, std::string_view ref_seq);

  Merge& merge_;
  SampleDAG sample_dag_;
  bool collapse_empty_fragment_edges_;
  ArbitraryInt applied_moves_count_;
  ConcurrentUnorderedMap<std::string, CompactGenome> sample_id_to_cg_map_;
  ReassignedStatesStorage reassigned_states_storage_ =
      AddMappedNodes(AddMATConversion(Storage{{}}));
  std::shared_mutex mat_mtx_;
  std::unique_ptr<MATStorage> sample_mat_storage_;
  std::mutex merge_mtx_;
  tbb::concurrent_vector<SPRType> batch_storage_;
  tbb::concurrent_vector<FragmentStorageFor<decltype(std::declval<SPRType>().View())>>
      batch_;
};

#include "larch/impl/spr/batching_callback_impl.hpp"

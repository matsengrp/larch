
template <typename CRTP, typename SampleDAG>
BatchingCallback<CRTP, SampleDAG>::BatchingCallback(Merge& merge, SampleDAG sample_dag)
    : merge_{merge}, sample_dag_{sample_dag}, collapse_empty_fragment_edges_{true} {}
template <typename CRTP, typename SampleDAG>
BatchingCallback<CRTP, SampleDAG>::BatchingCallback(Merge& merge, SampleDAG sample_dag,
                                                    bool collapse_empty_fragment_edges)
    : merge_{merge},
      sample_dag_{sample_dag},
      collapse_empty_fragment_edges_{collapse_empty_fragment_edges} {}

template <typename CRTP, typename SampleDAG>
bool BatchingCallback<CRTP, SampleDAG>::operator()(
    Profitable_Moves& move, int best_score_change,
    std::vector<Node_With_Major_Allele_Set_Change>&
        nodes_with_major_allele_set_change) {
  Assert(move.src != nullptr);
  Assert(move.dst != nullptr);
  auto& storage = [this]() -> SPRType& {
    std::shared_lock lock{mat_mtx_};
    Assert(sample_mat_storage_ != nullptr);
    return *batch_storage_.push_back(SPRStorage(sample_mat_storage_->View()));
  }();

  storage.View().GetRoot().Validate(true);

  if (storage.View().InitHypotheticalTree(move, nodes_with_major_allele_set_change)) {
    storage.View().GetRoot().Validate(true);
    auto fragment = collapse_empty_fragment_edges_
                        ? storage.View().MakeFragment()
                        : storage.View().MakeUncollapsedFragment();

    auto& impl = static_cast<CRTP&>(*this);
    std::pair<bool, bool> accepted =
        impl.OnMove(storage.View(), fragment, move, best_score_change,
                    nodes_with_major_allele_set_change);

    if (accepted.first) {
      applied_moves_count_++;
      batch_.push_back(std::move(fragment));
      /*
      if (batch_.size() > 2048) {
        std::unique_lock lock{merge_mtx_};
        if (batch_.size() > 2048) {
          merge_.AddDAGs(batch_);
          merge_.GetResult().GetRoot().Validate(true, true);
          batch_.clear();
          batch_storage_.clear();
        }
      }
      */
    }

    return accepted.second;

  } else {
    return false;
  }
}

template <typename CRTP, typename SampleDAG>
void BatchingCallback<CRTP, SampleDAG>::operator()(MAT::Tree& tree) {
  std::cout << "Larch-Usher callback Applying " << applied_moves_count_ << "\n"
            << std::flush;
  applied_moves_count_ = 0;
  reassigned_states_storage_ = AddMappedNodes(AddMATConversion(Storage{{}}));
  reassigned_states_storage_.View().BuildFromMAT(
      tree, merge_.GetResult().GetReferenceSequence());
  check_edge_mutations(reassigned_states_storage_.View().Const());
  reassigned_states_storage_.View().RecomputeCompactGenomes(true);
  reassigned_states_storage_.View().SampleIdsFromCG();
  {
    std::unique_lock lock{merge_mtx_};
    if (not batch_.empty()) {
      merge_.AddDAGs(batch_);
      batch_.clear();
      batch_storage_.clear();
    }
    merge_.AddDAGs(std::vector{reassigned_states_storage_.View()});
    merge_.GetResult().GetRoot().Validate(true, true);
    merge_.ComputeResultEdgeMutations();
  }
  {
    std::unique_lock lock{mat_mtx_};
    CreateMATStorage(tree, merge_.GetResult().GetReferenceSequence());
  }
  static_cast<CRTP&>(*this).OnRadius();
}

template <typename CRTP, typename SampleDAG>
void BatchingCallback<CRTP, SampleDAG>::OnReassignedStates(MAT::Tree& tree) {
  applied_moves_count_ = 0;
  reassigned_states_storage_.View().BuildFromMAT(
      tree, merge_.GetResult().GetReferenceSequence());
  check_edge_mutations(reassigned_states_storage_.View().Const());
  reassigned_states_storage_.View().RecomputeCompactGenomes(true);
  reassigned_states_storage_.View().SampleIdsFromCG();
  {
    std::unique_lock lock{merge_mtx_};
    merge_.AddDAGs(std::vector{reassigned_states_storage_.View()});
    merge_.GetResult().GetRoot().Validate(true, true);
    merge_.ComputeResultEdgeMutations();
  }
  {
    std::unique_lock lock{mat_mtx_};
    CreateMATStorage(tree, merge_.GetResult().GetReferenceSequence());
  }
}

template <typename CRTP, typename SampleDAG>
Merge& BatchingCallback<CRTP, SampleDAG>::GetMerge() {
  return merge_;
}

template <typename CRTP, typename SampleDAG>
ArbitraryInt BatchingCallback<CRTP, SampleDAG>::GetAppliedMovesCount() {
  return applied_moves_count_;
}

template <typename CRTP, typename SampleDAG>
auto BatchingCallback<CRTP, SampleDAG>::GetMappedStorage() {
  return reassigned_states_storage_.View();
}

template <typename CRTP, typename SampleDAG>
void BatchingCallback<CRTP, SampleDAG>::CreateMATStorage(MAT::Tree& tree,
                                                         std::string_view ref_seq) {
  sample_mat_storage_ = std::make_unique<MATStorage>(AddMATConversion(Storage{{}}));
  auto view = sample_mat_storage_->View();
  view.BuildFromMAT(tree, ref_seq);
  check_edge_mutations(view.Const());
  view.RecomputeCompactGenomes(true);
  view.SampleIdsFromCG();
}

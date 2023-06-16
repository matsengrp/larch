
template <typename CRTP>
bool BatchingCallback<CRTP>::operator()(Profitable_Moves& move, int best_score_change,
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
    auto fragment = storage.View().MakeFragment();
    batch_.push_back(std::move(fragment));

    if (batch_.size() > 2048) {
      std::unique_lock lock{merge_mtx_};
      if (batch_.size() > 2048) {
        merge_.AddDAGs(batch_);
        merge_.GetResult().GetRoot().Validate(true, true);
        batch_.clear();
        batch_storage_.clear();
      }
    }
    auto& impl = static_cast<CRTP&>(*this);
    return impl.OnMove(storage.View(), fragment, move, best_score_change,
                       nodes_with_major_allele_set_change);
  } else {
    return false;
  }
}

template <typename CRTP>
void BatchingCallback<CRTP>::operator()(MAT::Tree& tree) {
  auto storage = AddMATConversion(Storage{{}});
  storage.View().BuildFromMAT(tree, merge_.GetResult().GetReferenceSequence());
  storage.View().RecomputeCompactGenomes(true);
  {
    std::unique_lock lock{merge_mtx_};
    if (not batch_.empty()) {
      merge_.AddDAGs(batch_);
      batch_.clear();
      batch_storage_.clear();
    }
    merge_.AddDAGs(std::vector{storage.View()});
    merge_.GetResult().GetRoot().Validate(true, true);
    merge_.ComputeResultEdgeMutations();
  }
  {
    std::unique_lock lock{mat_mtx_};
    CreateMATStorage(tree, merge_.GetResult().GetReferenceSequence());
  }
}

template <typename CRTP>
void BatchingCallback<CRTP>::OnReassignedStates(MAT::Tree& tree) {
  reassigned_states_storage_.View().BuildFromMAT(
      tree, merge_.GetResult().GetReferenceSequence());
  check_edge_mutations(reassigned_states_storage_.View().Const());
  reassigned_states_storage_.View().RecomputeCompactGenomes(true);
  {
    std::unique_lock lock{merge_mtx_};
    merge_.AddDAGs(std::vector{reassigned_states_storage_.View()});
    merge_.GetResult().GetRoot().Validate(true, true);
    merge_.ComputeResultEdgeMutations();
  }
}

template <typename CRTP>
void BatchingCallback<CRTP>::CreateMATStorage(MAT::Tree& tree,
                                              std::string_view ref_seq) {
  sample_mat_storage_ = std::make_unique<MATStorage>(AddMATConversion(Storage{{}}));
  auto view = sample_mat_storage_->View();
  view.BuildFromMAT(tree, ref_seq);
  check_edge_mutations(view.Const());
  view.RecomputeCompactGenomes(true);
}

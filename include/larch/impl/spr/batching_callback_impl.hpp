
template <typename CRTP, typename SampleDAG>
BatchingCallback<CRTP, SampleDAG>::BatchingCallback(Merge& merge, SampleDAG sample_dag)
    : merge_{merge}, sample_dag_{sample_dag} {}

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
    return batch_storage_
        .Emplace(std::this_thread::get_id(), SPRStorage(sample_mat_storage_->View()))
        .first;
  }();

  storage.View().GetRoot().Validate(true);

  if (storage.View().InitHypotheticalTree(move, nodes_with_major_allele_set_change)) {
    storage.View().GetRoot().Validate(true);
    auto fragment = storage.View().MakeFragment();

    auto& impl = static_cast<CRTP&>(*this);
    std::pair<bool, bool> accepted =
        impl.OnMove(storage.View(), fragment, move, best_score_change,
                    nodes_with_major_allele_set_change);

    if (accepted.first) {
      if (batch_.Emplace(std::this_thread::get_id(), std::move(fragment)).second >
          std::thread::hardware_concurrency()) {
        // TODO view instead of vector
        std::vector<decltype(std::declval<FragmentStorage<SPRViewType>>().View())>
            batch_vec;
        batch_.Consume([&](auto rng) {
          for (auto& i : rng) {
            batch_vec.push_back(i.View());
          }
        });
        // std::lock_guard lock{merge_mtx_};
        // merge_.AddDAGs(batch_vec | ranges::views::all);
        // std::unique_lock lock{mat_mtx_};
        // batch_storage_.Clear();
      }
    }

    return accepted.second;

  } else {
    return false;
  }
}

template <typename CRTP, typename SampleDAG>
void BatchingCallback<CRTP, SampleDAG>::operator()(MAT::Tree& tree) {
  reassigned_states_storage_ = AddMappedNodes(AddMATConversion(Storage{{}}));
  reassigned_states_storage_.View().BuildFromMAT(
      tree, merge_.GetResult().GetReferenceSequence());
  check_edge_mutations(reassigned_states_storage_.View().Const());
  reassigned_states_storage_.View().RecomputeCompactGenomes(true);
  {
    std::lock_guard lock{merge_mtx_};
    batch_.Consume([&](auto rng) {
      if (not rng.empty()) {
        // TODO view instead of vector
        std::vector<decltype(std::declval<FragmentStorage<SPRViewType>>().View())>
            batch_vec;
        for (auto& i : rng) {
          batch_vec.push_back(i.View());
        }
        merge_.AddDAGs(batch_vec | ranges::views::all);
        // batch_storage_.Clear();
      }
    });
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
  reassigned_states_storage_.View().BuildFromMAT(
      tree, merge_.GetResult().GetReferenceSequence());
  check_edge_mutations(reassigned_states_storage_.View().Const());
  reassigned_states_storage_.View().RecomputeCompactGenomes(true);
  {
    std::lock_guard lock{merge_mtx_};
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
}

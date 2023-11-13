
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
    return *batch_storage_.push_back(AddSPRStorage(sample_mat_storage_->View()));
  }();

  storage.View().GetRoot().Validate(true);
  for (auto i : storage.View().Const().GetLeafs()) {
    Assert(i.HaveSampleId());
  }

  if (storage.View().InitHypotheticalTree(move, nodes_with_major_allele_set_change)) {
    storage.View().GetRoot().Validate(true);

    auto fragment_storage = collapse_empty_fragment_edges_
                                ? storage.View().MakeFragment()
                                : storage.View().MakeUncollapsedFragment();
    auto fragment = fragment_storage.View();

    auto& impl = static_cast<CRTP&>(*this);
    std::pair<bool, bool> accepted =
        impl.OnMove(storage.View(), fragment, move, best_score_change,
                    nodes_with_major_allele_set_change);
    if (accepted.first) {
      for (auto node : fragment.GetNodes()) {
        if (not(node.IsUA() or node.IsMoveNew())) {
          Assert(node.GetId().value != NoId);
          if (node.GetOld().IsLeaf()) {
            Assert(node.GetOld().HaveSampleId());
            Assert(not node.GetOld().GetSampleId().value().empty());
          }
        }
      }
      for (auto edge : fragment.GetEdges()) {
        Assert(edge.GetId().value != NoId);
        Assert(edge.GetChild().GetId().value != NoId);
        if (not edge.GetParent().IsUA()) {
          Assert(edge.GetParent().GetId().value != NoId);
        }
      }
      applied_moves_count_++;
      batch_.push_back(std::move(fragment_storage));
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
  reassigned_states_storage_ = std::make_unique<ReassignedStatesStorage>(
      AddMappedNodes(AddMATConversion(Storage::EmptyDefault())));
  auto reassigned_states = reassigned_states_storage_->View();
  reassigned_states.BuildFromMAT(tree, merge_.GetResult().GetReferenceSequence());
  check_edge_mutations(reassigned_states.Const());
  reassigned_states.RecomputeCompactGenomes(true);
  {
    std::unique_lock lock{merge_mtx_};
    if (not batch_.empty()) {
      merge_.AddDAGs(batch_ | Transform::ToView());
      batch_.clear();
      batch_storage_.clear();
    }
    merge_.AddDAGs(std::vector{reassigned_states});
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
  Assert(reassigned_states_storage_);
  auto reassigned_states = reassigned_states_storage_->View();
  reassigned_states.BuildFromMAT(tree, merge_.GetResult().GetReferenceSequence());
  reassigned_states.RecomputeCompactGenomes();
  check_edge_mutations(reassigned_states.Const());
  reassigned_states.RecomputeCompactGenomes(false);
  {
    std::unique_lock lock{merge_mtx_};
    merge_.AddDAGs(std::vector{reassigned_states});
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
  Assert(reassigned_states_storage_);
  return reassigned_states_storage_->View();
}

template <typename CRTP, typename SampleDAG>
void BatchingCallback<CRTP, SampleDAG>::CreateMATStorage(MAT::Tree& tree,
                                                         std::string_view ref_seq) {
  sample_mat_storage_ =
      std::make_unique<MATStorage>(AddMATConversion(Storage::EmptyDefault()));
  auto view = sample_mat_storage_->View();
  view.BuildFromMAT(tree, ref_seq);
  check_edge_mutations(view.Const());
  view.RecomputeCompactGenomes(true);
  view.SampleIdsFromCG();
}

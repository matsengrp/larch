#include "larch/benchmark.hpp"

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

  finally merge_large_batch([&]() {
    if (moves_batch_.size_approx() > 100) {
      auto all = moves_batch_.GatherAndClear([this](auto buckets) {
        std::vector<MoveStorage> result;
        result.reserve(moves_batch_.size_approx());
        for (auto&& bucket : buckets) {
          for (auto&& stored_move : bucket) {
            Assert(stored_move.fragment);
            // GetFullDAG(stored_move.fragment->View()).GetRoot().Validate(true, false);
            result.push_back(std::move(stored_move));
          }
        }
        return result;
      });
      if (not all.empty()) {
        std::unique_lock lock{merge_mtx_};
        merge_.AddDAGs(
            all | ranges::views::transform([](auto& i) { return i.fragment->View(); }));
        // merge_.GetResult().GetRoot().Validate(true, true);
      }
    }
  });

  return moves_batch_.AddElement(
      [this, &move, best_score_change,
       &nodes_with_major_allele_set_change](auto& bucket) -> bool {
        std::shared_lock lock{mat_mtx_};

        /* CONDENSING CODE: probably want to change this to an uncondensed storage, once
         * it's implemented*/
        // TODO old
        Assert(sample_mat_storage_ != nullptr);
        bucket.push_back(MoveStorage{
            std::make_unique<SPRType>(AddSPRStorage(sample_mat_storage_->View())),
            nullptr});
        // TODO new
        // Assert(sample_matview_storage_ != nullptr);
        // bucket.push_back(MoveStorage{
        //     std::make_unique<NewSPRType>(AddSPRStorage(sample_matview_storage_->View())),
        //     nullptr});

        auto& storage = bucket.back();
    // storage.spr->View().GetRoot().Validate(true);
#ifndef NDEBUG
        for (auto i : storage.spr->View().Const().GetLeafs()) {
          Assert(i.HaveSampleId());
        }
#endif
        if (storage.spr->View().InitHypotheticalTree(
                move, nodes_with_major_allele_set_change)) {
          // storage.spr->View().GetRoot().Validate(true);
          storage.fragment = std::make_unique<FragmentType>(
              collapse_empty_fragment_edges_
                  ? storage.spr->View().MakeFragment()
                  : storage.spr->View().MakeUncollapsedFragment());
          auto fragment = storage.fragment->View();
          // GetFullDAG(fragment).GetRoot().Validate(true, false);
          auto& impl = static_cast<CRTP&>(*this);
          std::pair<bool, bool> accepted =
              impl.OnMove(storage.spr->View(), fragment, move, best_score_change,
                          nodes_with_major_allele_set_change);
          if (accepted.first) {
#ifndef NDEBUG
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
#endif
            applied_moves_count_.fetch_add(1);
          }

          return accepted.second;

        } else {
          bucket.pop_back();
          return false;
        }
      });
}

template <typename CRTP, typename SampleDAG>
void BatchingCallback<CRTP, SampleDAG>::operator()(MAT::Tree& tree) {
  std::cout << "Larch-Usher callback Applying " << applied_moves_count_.load() << "\n"
            << std::flush;
  applied_moves_count_.store(0);

  // TODO old
  reassigned_states_storage_ = std::make_unique<ReassignedStatesStorage>(
      AddMappedNodes(AddMATConversion(Storage::EmptyDefault())));
  // TODO new
  // new_reassigned_states_storage_ = std::make_unique<NewReassignedStatesStorage>(
  //     AddMappedNodes(AddMATConversion(NewStorage::EmptyDefault())));

  auto reassigned_states = reassigned_states_storage_->View();
  reassigned_states.BuildFromMAT(tree, merge_.GetResult().GetReferenceSequence());
  check_edge_mutations(reassigned_states.Const());
  reassigned_states.RecomputeCompactGenomes(true);
  {
    auto all = moves_batch_.GatherAndClear([this](auto buckets) {
      std::vector<MoveStorage> result;
      result.reserve(moves_batch_.size_approx());
      for (auto&& bucket : buckets) {
        for (auto&& stored_move : bucket) {
          Assert(stored_move.fragment);
          // GetFullDAG(stored_move.fragment->View()).GetRoot().Validate(true, false);
          result.push_back(std::move(stored_move));
        }
      }
      return result;
    });
    std::unique_lock lock{merge_mtx_};
    if (not all.empty()) {
      merge_.AddDAGs(
          all | ranges::views::transform([](auto& i) { return i.fragment->View(); }));
    }
    merge_.AddDAGs(std::vector{reassigned_states});
    // merge_.GetResult().GetRoot().Validate(true, true);
  }
  {
    std::unique_lock lock{mat_mtx_};
    // TODO old
    CreateMATStorage(tree, merge_.GetResult().GetReferenceSequence());
    // TODO new
    CreateMATViewStorage(tree, merge_.GetResult().GetReferenceSequence());
  }
  static_cast<CRTP&>(*this).OnRadius();
}

template <typename CRTP, typename SampleDAG>
void BatchingCallback<CRTP, SampleDAG>::OnReassignedStates(MAT::Tree& tree) {
  applied_moves_count_.store(0);
  Assert(reassigned_states_storage_);
  auto reassigned_states = reassigned_states_storage_->View();
  reassigned_states.BuildFromMAT(tree, merge_.GetResult().GetReferenceSequence());
  reassigned_states.RecomputeCompactGenomes();
  check_edge_mutations(reassigned_states.Const());
  reassigned_states.RecomputeCompactGenomes(false);
  {
    std::unique_lock lock{merge_mtx_};
    merge_.AddDAGs(std::vector{reassigned_states});
    // merge_.GetResult().GetRoot().Validate(true, true);
  }
  {
    std::unique_lock lock{mat_mtx_};
    // TODO old
    CreateMATStorage(tree, merge_.GetResult().GetReferenceSequence());
    // TODO new
    CreateMATViewStorage(tree, merge_.GetResult().GetReferenceSequence());
  }
}

template <typename CRTP, typename SampleDAG>
Merge& BatchingCallback<CRTP, SampleDAG>::GetMerge() {
  return merge_;
}

template <typename CRTP, typename SampleDAG>
size_t BatchingCallback<CRTP, SampleDAG>::GetAppliedMovesCount() {
  return applied_moves_count_.load();
}

template <typename CRTP, typename SampleDAG>
auto BatchingCallback<CRTP, SampleDAG>::GetMappedStorage() {
  Assert(reassigned_states_storage_);
  return reassigned_states_storage_->View();
}

// TODO old
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

// TODO new
template <typename CRTP, typename SampleDAG>
void BatchingCallback<CRTP, SampleDAG>::CreateMATViewStorage(MAT::Tree& tree,
                                                             std::string_view ref_seq) {
  std::ignore = tree;
  std::ignore = ref_seq;
  sample_matview_storage_ = std::make_unique<MATStorageImpl>();
  sample_matview_storage_->GetCondensed().View().SetMAT(std::addressof(tree));
  auto storage =
      NewStorage::Consume(std::move(sample_matview_storage_->GetCondensed()));
  auto view = storage.View();
  view.SetReferenceSequence(ref_seq);
  view.BuildRootAndLeafs();
  view.RecomputeCompactGenomes(true);
}

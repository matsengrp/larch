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
      Benchmark bench;
      bench.start();
      auto all = moves_batch_.GatherAndClear([this](auto buckets) {
        std::vector<MoveStorage> result;
        result.reserve(moves_batch_.size_approx());
        for (auto&& bucket : buckets) {
          for (auto&& stored_move : bucket) {
            if (stored_move.fragment) {
              result.push_back(std::move(stored_move));
            }
          }
        }
        return result;
      });
      std::cout << "Gather " << all.size() << " in " << bench.lapMs() << "ms\n";
      std::unique_lock lock{merge_mtx_};
      bench.lapMs();
      merge_.AddDAGs(
          all | ranges::views::transform([](auto& i) { return i.fragment->View(); }));
      std::cout << "  Merged in " << bench.lapMs() << "ms\n";
      // merge_.GetResult().GetRoot().Validate(true, true);
    }
  });

  return moves_batch_.AddElement(
      [this, &move, best_score_change,
       &nodes_with_major_allele_set_change](auto& bucket) -> bool {
        Benchmark bench;
        std::shared_lock lock{mat_mtx_};
        Assert(sample_mat_storage_ != nullptr);
        bench.start();
        bucket.push_back(MoveStorage{
            std::make_unique<SPRType>(AddSPRStorage(sample_mat_storage_->View())),
            nullptr});
        // std::cout << "Add spr: " << bench.lapMs() << "ms\n";
        auto& storage = bucket.back();
    // storage.spr->View().GetRoot().Validate(true);
#ifndef NDEBUG
        for (auto i : storage.spr->View().Const().GetLeafs()) {
          Assert(i.HaveSampleId());
        }
#endif
        bench.lapMs();
        if (storage.spr->View().InitHypotheticalTree(
                move, nodes_with_major_allele_set_change)) {
          // storage.spr->View().GetRoot().Validate(true);
          // std::cout << "  Init: " << bench.lapMs() << "ms\n";
          storage.fragment = std::make_unique<FragmentType>(
              collapse_empty_fragment_edges_
                  ? storage.spr->View().MakeFragment()
                  : storage.spr->View().MakeUncollapsedFragment());
          // std::cout << "  Fragment: " << bench.lapMs() << "ms\n";
          auto fragment = storage.fragment->View();

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
          return false;
        }
      });
}

template <typename CRTP, typename SampleDAG>
void BatchingCallback<CRTP, SampleDAG>::operator()(MAT::Tree& tree) {
  std::cout << "Larch-Usher callback Applying " << applied_moves_count_.load() << "\n"
            << std::flush;
  applied_moves_count_.store(0);
  reassigned_states_storage_ = std::make_unique<ReassignedStatesStorage>(
      AddMappedNodes(AddMATConversion(Storage::EmptyDefault())));
  auto reassigned_states = reassigned_states_storage_->View();
  reassigned_states.BuildFromMAT(tree, merge_.GetResult().GetReferenceSequence());
  check_edge_mutations(reassigned_states.Const());
  reassigned_states.RecomputeCompactGenomes(true);
  {
    std::unique_lock lock{merge_mtx_};
    if (moves_batch_.size_approx() > 0) {
      auto all = moves_batch_.GatherAndClear([this](auto buckets) {
        std::vector<MoveStorage> result;
        result.reserve(moves_batch_.size_approx());
        for (auto&& bucket : buckets) {
          for (auto&& stored_move : bucket) {
            if (stored_move.fragment) {
              result.push_back(std::move(stored_move));
            }
          }
        }
        return result;
      });
      merge_.AddDAGs(
          all | ranges::views::transform([](auto& i) { return i.fragment->View(); }));
    }
    merge_.AddDAGs(std::vector{reassigned_states});
    // merge_.GetResult().GetRoot().Validate(true, true);
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
size_t BatchingCallback<CRTP, SampleDAG>::GetAppliedMovesCount() {
  return applied_moves_count_.load();
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

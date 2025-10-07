#include "larch/benchmark.hpp"

template <typename CRTP>
BatchingCallback<CRTP>::BatchingCallback(Merge& merge)
    : merge_{merge}, collapse_empty_fragment_edges_{true} {}
template <typename CRTP>
BatchingCallback<CRTP>::BatchingCallback(Merge& merge,
                                         bool collapse_empty_fragment_edges)
    : merge_{merge}, collapse_empty_fragment_edges_{collapse_empty_fragment_edges} {}

template <typename CRTP>
bool BatchingCallback<CRTP>::operator()(Profitable_Moves& move, int best_score_change,
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

        Assert(sample_mat_storage_ != nullptr);
        bucket.push_back(MoveStorage{
            std::make_unique<SPRType>(AddSPRStorage(sample_mat_storage_->View())),
            nullptr});

        auto& storage = bucket.back();
    // MADAGToDOT(storage.spr->View(), std::cout);
    // storage.spr->View().GetRoot().Validate(true);
#ifdef KEEP_ASSERTS
        for (auto i : storage.spr->View().Const().GetLeafs()) {
          Assert(i.HaveSampleId());
        }
        auto old_src_parent = storage.spr->View().GetNodeFromMAT(move.src->parent);
#endif
        if (storage.spr->View().InitHypotheticalTree(
                move, nodes_with_major_allele_set_change)) {
          // storage.spr->View().GetRoot().Validate(true);
          storage.fragment = std::make_unique<FragmentType>(
              collapse_empty_fragment_edges_
                  ? storage.spr->View().MakeFragment()
                  : storage.spr->View().MakeUncollapsedFragment());
#ifdef KEEP_ASSERTS
          auto new_old_src_parent = storage.spr->View().GetOldSourceParent();
          Assert(old_src_parent.GetId() == new_old_src_parent.GetId());
#endif
          auto fragment = storage.fragment->View();
          // GetFullDAG(fragment).GetRoot().Validate(true, false);
          auto& impl = static_cast<CRTP&>(*this);
          std::pair<bool, bool> accepted =
              impl.OnMove(storage.spr->View(), fragment, move, best_score_change,
                          nodes_with_major_allele_set_change);

          if (accepted.first) {
#ifdef KEEP_ASSERTS
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
          } else {
            bucket.pop_back();
          }

          return accepted.second;

        } else {
          bucket.pop_back();
          return false;
        }
      });
}

template <typename CRTP>
void BatchingCallback<CRTP>::operator()(MAT::Tree& tree) {
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
#ifdef KEEP_ASSERTS
      auto orig_num_leafs = merge_.GetResult().GetLeafsCount();
#endif
      merge_.AddDAGs(
          all | ranges::views::transform([](auto& i) { return i.fragment->View(); }));

      //// FOR DEBUGGING: alternative "serialized" merging of fragments to check them
      // for (auto& i : all) {
      //   auto frag = i.fragment->View();
      //   MADAGToDOT(frag, std::cout);
      //   MADAGToDOT(i.spr->View(), std::cout);
      //   std::vector dags = {frag};
      //   merge_.AddDAGs(dags);
      //   MADAGToDOT(merge_.GetResult(), std::cout);
      // Assert(merge_.GetResult().GetLeafsCount() == orig_num_leafs);
      // }

#ifdef KEEP_ASSERTS
      Assert(merge_.GetResult().GetLeafsCount() == orig_num_leafs);
#endif
    }
    merge_.AddDAGs(std::vector{reassigned_states});
    // merge_.GetResult().GetRoot().Validate(true, true);
  }
  {
    std::unique_lock lock{mat_mtx_};
#if USE_MAT_VIEW
    CreateMATViewStorage(tree, merge_.GetResult().GetReferenceSequence());
#else
    CreateMATStorage(tree, merge_.GetResult().GetReferenceSequence());
#endif
  }
  static_cast<CRTP&>(*this).OnRadius();
}

template <typename CRTP>
void BatchingCallback<CRTP>::OnReassignedStates(MAT::Tree& tree) {
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
#if USE_MAT_VIEW
    CreateMATViewStorage(tree, merge_.GetResult().GetReferenceSequence());
#else
    CreateMATStorage(tree, merge_.GetResult().GetReferenceSequence());
#endif
  }
}

template <typename CRTP>
Merge& BatchingCallback<CRTP>::GetMerge() {
  return merge_;
}

template <typename CRTP>
size_t BatchingCallback<CRTP>::GetAppliedMovesCount() {
  return applied_moves_count_.load();
}

template <typename CRTP>
auto BatchingCallback<CRTP>::GetMappedStorage() {
  Assert(reassigned_states_storage_);
  return reassigned_states_storage_->View();
}

#if USE_MAT_VIEW
template <typename CRTP>
void BatchingCallback<CRTP>::CreateMATViewStorage(MAT::Tree& tree,
                                                  std::string_view ref_seq) {
  SetSample(tree, std::string{ref_seq});
  UncondensedMATViewStorage mv_storage;
  mv_storage.View().SetMAT(std::addressof(sample_mat_tree_));
  // mv_storage.View().BuildRootAndLeafs();
  // MADAGToDOT(mv_storage.View(), std::cout);
  // mv_storage.View().GetRoot().Validate(true, false);

  sample_mat_storage_ = std::make_unique<MATStorage>(std::move(mv_storage));
  UncondensedMergeDAGStorage& storage = *sample_mat_storage_;
  auto view = storage.View();
  view.SetReferenceSequence(sample_refseq_);
  view.BuildRootAndLeafs();

  // alternative sampleId placement(please change if there's a better way!!)
  for (auto node : view.GetNodes()) {
    if (node.IsLeaf()) {
      std::string sample_id = "";
      if (node.GetMATNode() == nullptr) {
        auto& node_storage =
            node.template GetFeatureExtraStorage<MATNodeStorage>().get();
        sample_id = node_storage.node_id_to_sampleid_map_.at(node.GetId());
      } else {
        sample_id = sample_mat_tree_.get_node_name(node.GetMATNode()->node_id);
      }
      Assert(not sample_id.empty());
      auto id_iter = view.template AsFeature<Deduplicate<SampleId>>().AddDeduplicated(
          SampleId::Make(std::move(sample_id)));
      node = id_iter.first;
    }
  }
  /*
    for (auto node : view.GetNodes()) {
      if (node.GetMATNode() == nullptr) {
        continue;
      }
      std::string sample_id = tree.condensed_nodes.at(node.GetMATNode()->node_id).at(0);
      Assert(not sample_id.empty());
      auto id_iter = view.template AsFeature<Deduplicate<SampleId>>().AddDeduplicated(
          SampleId{std::move(sample_id)});
      node = id_iter.first;
    }
  */
  view.RecomputeCompactGenomes(true);
}
#else
template <typename CRTP>
void BatchingCallback<CRTP>::CreateMATStorage(MAT::Tree& tree,
                                              std::string_view ref_seq) {
  sample_mat_storage_ =
      std::make_unique<MATStorage>(AddMATConversion(Storage::EmptyDefault()));
  auto view = sample_mat_storage_->View();
  view.BuildFromMAT(tree, ref_seq);
  check_edge_mutations(view.Const());
  view.RecomputeCompactGenomes(true);
  view.SampleIdsFromCG();
}
#endif

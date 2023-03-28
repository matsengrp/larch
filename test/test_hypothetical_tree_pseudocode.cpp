template <typename Fragment, typename SampleDAG, typename parsimonyScoreType>
parsimonyScoreType ComputeParsimonyScoreChangeForFragment(Fragment fragment,
                                                          SampleDAG sample) {
  // calculate the parsimony score for the given fragment using sampleDAG to read
  // compact genomes

  parsimonyScoreType orig_parsimony = 0;
  parsimonyScoreType new_parsimony = 0;
  // traverse the nodes in the fragment and find the corresponding parent/child pair in
  // sampleDAG. Since the fragment has an altered topology from sampleDAG, we retrieve
  // the corresponding sampleDAG's node. We then compute the CG hamming distance for the
  // sampleDAG's edge above that node, as well as the CG hamming distance for the
  // fragment's edge above that node. we do not have to worry that

  for (auto fragment_child_node : fragment) {
    auto sample_child_node = sample.Get(fragment_child_node);
    auto sample_parent_node = sample_child_node.GetSingleParent().GetParent();
    if (!sample_parent_node.IsRoot()) {
      auto fragment_parent_node = fragment_child_node.GetSingleParent().GetParent();
      orig_parsimony += sample_child_node.GetCompactGenome()
                            .DifferingSites(sample_parent_node.GetCompactGenome())
                            .size();
      new_parsimony += fragment_child_node.GetCompactGenome()
                           .DifferingSites(fragment_parent_node.GetCompactGenome())
                           .size();
    }
  }

  return new_parsimony - orig_parsimony;
}

template <typename SampleDAG, typename parsimonyScoreType>
struct Single_Move_Callback_With_Hypothetical_Tree : public Move_Found_Callback {
  Single_Move_Callback_With_Hypothetical_Tree(Merge<MADAG>& merge, SampleDAG sample)
      : merge_{merge},
        sample_{sample},
        approved_a_move_{false},
        computed_score_change_{0.0},
        matOptimize_reported_score_change{0.0} {}

  bool operator()(Profitable_Moves& move, int /*best_score_change*/,
                  std::vector<Node_With_Major_Allele_Set_Change>&
                      nodes_with_major_allele_set_change) override {
    if (!approved_a_move_) {
      // apply move to merge object.

      auto spr_storage = SPRStorage(sample_);
      auto spr = spr_storage.View();

      // ** create hypothetical tree
      spr.InitHypotheticalTree(move, nodes_with_major_allele_set_change);

      // ** build fragment
      auto spr_fragment = spr.GetFragment();

      // calculate parsimony score change for this fragment and save, along with
      // matOptimize's computed score change
      computed_score_change_ =
          ComputeParsimonyScoreChangeForFragment(spr_fragment, sample_);
      matOptimize_reported_score_change_ = move.score_change;

      // set flag so we don't approve any more moves
      approved_a_move_ = true;

      // ** merge fragment into merge
      // TODO merge_.AddDAG(spr_fragment);

      // return true so we do apply this move.
      return true;

    } else {
      return !approved_a_move_;
    }
  }
  Merge<MADAG>& merge_;
  SampleDAG sample_;
  bool approved_a_move_;
  parsimonyScoreType computed_score_change_;
  parsimonyScoreType matOptimize_reported_score_change_;
};

static void test_optimizing_with_hypothetical_tree(
    const MADAGStorage& tree_shaped_dag) {
  // this test takes a tree and uses matOptimize to apply a single move.

  Merge<MADAG> dag_altered_in_callback{tree_shaped_dag.View().GetReferenceSequence()};
  dag_altered_in_callback.AddDAGs({tree_shaped_dag.View()});

  // sample tree
  SubtreeWeight<ParsimonyScore, MergeDAG> weight{dag_altered_in_callback.GetResult()};
  auto sample = AddMATConversion(weight.SampleTree({}));
  check_edge_mutations(sample.View());
  MAT::Tree mat;
  sample.View().BuildMAT(mat);

  // create a callback that only allows one move
  Single_Move_Callback_With_Hypothetical_Tree single_move_callback{
      dag_altered_in_callback, sample.View()};

  // optimize tree with matOptimize using a callback that only applies a single move
  auto optimized_tree =
      optimize_dag_direct(sample.View(), single_move_callback, [](MAT::Tree) {});

  Merge<MADAG> two_tree_dag{tree_shaped_dag.View().GetReferenceSequence()};
  two_tree_dag.AddDAG(sample.View());
  two_tree_dag.AddDAG(optimized_tree.first.View());

  // check topologies of the two DAGs match in feature count
  Assert(two_tree_dag.GetResult().GetNodesCount() ==
         dag_altered_in_callback.GetResult().GetNodesCount());
  Assert(two_tree_dag.GetResult().GetEdgesCount() ==
         dag_altered_in_callback.GetResult().GetEdgesCount());

  Assert(single_move_callback.computed_score_change ==
         single_move_callback.matOptimize_reported_score_change_);
}

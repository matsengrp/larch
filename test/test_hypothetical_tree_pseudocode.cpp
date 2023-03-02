struct Single_Move_Callback_With_Hypothetical_Tree : public Move_Found_Callback {
  Single_Move_Callback(const Merge<MADAG>& merge, SampleDAG sample)
      : merge_{merge},
        sample_{sample},
        approved_a_move_{false} {}

  bool operator()(Profitable_Moves& move, int best_score_change,
                  std::vector<Node_With_Major_Allele_Set_Change>& nodes_with_major_allele_set_change) override {
    if (!approved_a_move) {
      // apply move to merge object.

      auto spr_storage = SPRStorage(sample_);
      auto spr = spr_storage.View();

      // ** create hypothetical tree
      spr.InitHypotheticalTree(move, nodes_with_major_allele_set_change);
      spr.ApplyMove(move.src->node_id, move.dst->node_id);

      // ** build fragment
      auto spr_fragment = spr.GetFragment();

      // set flag so we don't approve any more moves
      approved_a_move_ = true;

      // ** merge fragment into merge
      merge.AddDAGs(spr_fragment);

      // return true so we do apply this move.
      return true;

    } else {
      return !approved_a_move_;
    }
  }
  Merge merge_;
  DAG sample_dag_;
  bool approved_a_move_;
};

static void test_optimizing_with_hypothetical_tree(std::string_view single_tree_path){
  // this test takes a tree and uses matOptimize to apply a single move.

  // load a tree shaped DAG
  MADAGStorage tree_shaped_dag = LoadDAGFromProtobuf(single_tree_path)
  tree_shaped_dag.RecomputeCompactGenomes();

  Merge<MADAG> dag_altered_in_callback{tree_shaped_dag.GetReferenceSequence()};
  dag_altered_in_callback.AddDAGs({tree_shaped_dag});

  // sample tree
  SubtreeWeight<ParsimonyScore, MergeDAG> weight{dag_altered_in_callback.GetResult()};
  auto sample = weight.SampleTree({}).first;
  check_edge_mutations(sample.View());

  // create a callback that only allows one move
  Single_Move_Callback_With_Hypothetical_Tree single_move_callback{dag_altered_in_callback, sample.View()};

  // optimize tree with matOptimize using a callback that only applies a single move
  auto optimized_tree = optimize_dag_direct(sample.View(), single_move_callback, [](MAT::Tree) {});

  Merge<MADAG> two_tree_dag{tree_shaped_dag.GetReferenceSequence()};
  two_tree_dag.AddDAGs(sample, optimized_tree);

  // check topologies of the two DAGs match in feature count
  Assert(two_tree_dag.GetNodesCount() == dag_altered_in_callback.GetNodesCount());
  Assert(two_tree_dag.GetEdges() == dag_altered_in_callback.GetEdgesCount());
}


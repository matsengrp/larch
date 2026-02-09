#include "larch/subtree/sankoff.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/subtree/subtree_weight.hpp"

#include "sample_dag.hpp"
#include "test_common.hpp"

// Test that BaseToIndex and IndexToBase are inverses
static void test_base_conversion() {
  TestAssert(BaseToIndex('A') == 0);
  TestAssert(BaseToIndex('C') == 1);
  TestAssert(BaseToIndex('G') == 2);
  TestAssert(BaseToIndex('T') == 3);

  TestAssert(IndexToBase(0) == 'A');
  TestAssert(IndexToBase(1) == 'C');
  TestAssert(IndexToBase(2) == 'G');
  TestAssert(IndexToBase(3) == 'T');

  // Round-trip
  for (char base : {'A', 'C', 'G', 'T'}) {
    TestAssert(IndexToBase(BaseToIndex(base)) == base);
  }
}

[[maybe_unused]] static const auto test_base_conversion_registered =
    add_test({test_base_conversion, "Sankoff: Base conversion", {"sankoff"}});

// Test GetCompatibleIndices for ambiguous and unambiguous bases
static void test_compatible_indices() {
  // Unambiguous bases
  auto a_indices = GetCompatibleIndices(MutationBase{'A'});
  TestAssert(a_indices.size() == 1);
  TestAssert(a_indices[0] == 0);

  auto c_indices = GetCompatibleIndices(MutationBase{'C'});
  TestAssert(c_indices.size() == 1);
  TestAssert(c_indices[0] == 1);

  // Ambiguous base 'N' should return all 4
  auto n_indices = GetCompatibleIndices(MutationBase{'N'});
  TestAssert(n_indices.size() == 4);
}

[[maybe_unused]] static const auto test_compatible_indices_registered =
    add_test({test_compatible_indices, "Sankoff: Compatible indices", {"sankoff"}});

// Test uniform cost matrix creation
static void test_uniform_cost_matrix() {
  SankoffScorer<MADAG> scorer(make_sample_dag().View());

  // Access the DP table to check cost matrices
  const auto& dp_table = scorer.GetDPTable();

  // All sites should have uniform cost matrices
  for (const auto& matrix : dp_table.cost_matrices) {
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        if (i == j) {
          TestAssert(matrix[i][j] == 0.0);
        } else {
          TestAssert(matrix[i][j] == 1.0);
        }
      }
    }
  }
}

[[maybe_unused]] static const auto test_uniform_cost_matrix_registered =
    add_test({test_uniform_cost_matrix, "Sankoff: Uniform cost matrix", {"sankoff"}});

// Test variable site detection
static void test_variable_site_detection() {
  auto dag_storage = make_sample_dag();
  auto dag = dag_storage.View();

  SankoffScorer<decltype(dag)> scorer(dag);

  // sample_dag has mutations at positions 1, 2, 3 (reference is "GAA")
  TestAssert(scorer.GetNumVariableSites() == 3);

  const auto& dp_table = scorer.GetDPTable();
  // Check all expected positions are in the map
  TestAssert(dp_table.site_to_index.count(1) > 0);
  TestAssert(dp_table.site_to_index.count(2) > 0);
  TestAssert(dp_table.site_to_index.count(3) > 0);
}

[[maybe_unused]] static const auto test_variable_site_detection_registered = add_test(
    {test_variable_site_detection, "Sankoff: Variable site detection", {"sankoff"}});

// Test basic scoring on sample DAG with uniform costs
static void test_basic_uniform_scoring() {
  auto dag_storage = make_sample_dag();
  auto dag = dag_storage.View();

  SankoffScorer<decltype(dag)> scorer(dag);

  // Compute score below root
  auto root = dag.GetRoot();
  double score = scorer.ComputeScoreBelow(root);

  // Sankoff with uniform costs finds the optimal reconstruction, so its score
  // should be at most the ParsimonyScore (which sums existing edge annotations).
  MADAG const_dag = dag;
  SubtreeWeight<ParsimonyScore, MADAG> weight{const_dag};
  ParsimonyScore::Weight parsimony =
      weight.ComputeWeightBelow(const_dag.GetRoot(), {});
  TestAssert(score <= static_cast<double>(parsimony));

  // Verify cached score matches
  TestAssert(scorer.GetTotalScore() == score);
}

[[maybe_unused]] static const auto test_basic_uniform_scoring_registered = add_test(
    {test_basic_uniform_scoring, "Sankoff: Basic uniform scoring", {"sankoff"}});

// Test ancestral reconstruction
static void test_ancestral_reconstruction() {
  auto dag_storage = make_sample_dag();
  auto dag = dag_storage.View();

  SankoffScorer<decltype(dag)> scorer(dag);

  auto root = dag.GetRoot();
  scorer.ComputeScoreBelow(root);
  scorer.ReconstructAncestralSequences(root);

  // Verify we can get reconstructed genomes for all nodes
  for (auto node : dag.GetNodes()) {
    auto genome = scorer.GetReconstructedGenome(node);
    // The genome should be valid (no exception thrown)
    // Can't easily verify content without knowing expected values
  }

  // Verify reconstructed bases are valid (0-3)
  size_t num_sites = scorer.GetNumVariableSites();
  for (auto node : dag.GetNodes()) {
    for (size_t site_idx = 0; site_idx < num_sites; ++site_idx) {
      char base = scorer.GetReconstructedBase(node.GetId(), site_idx);
      TestAssert(base == 'A' || base == 'C' || base == 'G' || base == 'T');
    }
  }
}

[[maybe_unused]] static const auto test_ancestral_reconstruction_registered = add_test(
    {test_ancestral_reconstruction, "Sankoff: Ancestral reconstruction", {"sankoff"}});

// Test that leaf nodes get their observed bases
static void test_leaf_base_assignment() {
  auto dag_storage = make_unambiguous_sample_dag();
  auto dag = dag_storage.View();

  SankoffScorer<decltype(dag)> scorer(dag);

  auto root = dag.GetRoot();
  scorer.ComputeScoreBelow(root);
  scorer.ReconstructAncestralSequences(root);

  // For leaf nodes, the reconstructed genome should match the original
  // CompactGenome at variable sites
  const std::string& ref_seq = dag.GetReferenceSequence();
  const auto& dp_table = scorer.GetDPTable();

  for (auto node : dag.GetNodes()) {
    if (node.IsLeaf()) {
      auto reconstructed = scorer.GetReconstructedGenome(node);
      const auto& original = node.GetCompactGenome();

      // For each variable site, check if reconstructed matches original
      for (size_t site_idx = 0; site_idx < scorer.GetNumVariableSites(); ++site_idx) {
        MutationPosition pos = dp_table.variable_sites[site_idx];
        char recon_base = scorer.GetReconstructedBase(node.GetId(), site_idx);

        // Get original base
        auto maybe_orig = original[pos];
        char orig_base;
        if (maybe_orig.has_value()) {
          orig_base = maybe_orig->ToChar();
        } else {
          orig_base = ref_seq.at(pos.value - 1);
        }

        TestAssert(recon_base == orig_base);
      }
    }
  }
}

[[maybe_unused]] static const auto test_leaf_base_assignment_registered =
    add_test({test_leaf_base_assignment, "Sankoff: Leaf base assignment", {"sankoff"}});

// Test with custom cost matrix (transition/transversion weighting)
static void test_transition_transversion_costs() {
  auto dag_storage = make_sample_dag();
  auto dag = dag_storage.View();

  // Create cost matrix where transitions (A<->G, C<->T) cost 1
  // and transversions cost 2
  SiteCostMatrix ts_tv_matrix;
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      if (i == j) {
        ts_tv_matrix[i][j] = 0.0;
      } else if ((i == 0 && j == 2) || (i == 2 && j == 0) ||  // A<->G
                 (i == 1 && j == 3) || (i == 3 && j == 1)) {  // C<->T
        ts_tv_matrix[i][j] = 1.0;                             // Transition
      } else {
        ts_tv_matrix[i][j] = 2.0;  // Transversion
      }
    }
  }

  SankoffScorer<decltype(dag)> scorer_uniform(dag);
  SankoffScorer<decltype(dag)> scorer_ts_tv(dag, ts_tv_matrix);

  auto root = dag.GetRoot();
  double score_uniform = scorer_uniform.ComputeScoreBelow(root);
  double score_ts_tv = scorer_ts_tv.ComputeScoreBelow(root);

  // With ts/tv weighting, score should be >= uniform score
  // (transversions cost more, transitions cost same)
  TestAssert(score_ts_tv >= score_uniform);
}

[[maybe_unused]] static const auto test_transition_transversion_costs_registered =
    add_test({test_transition_transversion_costs,
              "Sankoff: Transition/transversion costs",
              {"sankoff"}});

// Test with site-specific cost matrices
static void test_site_specific_costs() {
  auto dag_storage = make_sample_dag();
  auto dag = dag_storage.View();

  // Create different cost matrices for different sites
  std::unordered_map<size_t, SiteCostMatrix> site_matrices;

  // Site 1: uniform costs
  SiteCostMatrix uniform;
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      uniform[i][j] = (i == j) ? 0.0 : 1.0;
    }
  }
  site_matrices[1] = uniform;

  // Site 2: higher costs
  SiteCostMatrix high_cost;
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      high_cost[i][j] = (i == j) ? 0.0 : 5.0;
    }
  }
  site_matrices[2] = high_cost;

  // Site 3 uses default (uniform)

  SankoffScorer<decltype(dag)> scorer(dag, site_matrices);

  auto root = dag.GetRoot();
  double score = scorer.ComputeScoreBelow(root);

  // Score should be positive
  TestAssert(score > 0.0);

  // Verify the cost matrices were set correctly
  const auto& dp_table = scorer.GetDPTable();
  size_t site1_idx = dp_table.site_to_index.at(1);
  size_t site2_idx = dp_table.site_to_index.at(2);

  // Site 1 should have cost 1 for non-matching
  TestAssert(dp_table.cost_matrices[site1_idx][0][1] == 1.0);

  // Site 2 should have cost 5 for non-matching
  TestAssert(dp_table.cost_matrices[site2_idx][0][1] == 5.0);
}

[[maybe_unused]] static const auto test_site_specific_costs_registered =
    add_test({test_site_specific_costs, "Sankoff: Site-specific costs", {"sankoff"}});

// Test scoring with ambiguous bases at leaves
static void test_ambiguous_leaf_bases() {
  auto dag_storage = make_ambiguous_sample_dag();
  auto dag = dag_storage.View();

  SankoffScorer<decltype(dag)> scorer(dag);

  auto root = dag.GetRoot();
  double score = scorer.ComputeScoreBelow(root);

  // Should complete without error and produce a valid score
  TestAssert(score >= 0.0);

  // Reconstruct and verify all bases are unambiguous
  scorer.ReconstructAncestralSequences(root);

  size_t num_sites = scorer.GetNumVariableSites();
  for (auto node : dag.GetNodes()) {
    for (size_t site_idx = 0; site_idx < num_sites; ++site_idx) {
      char base = scorer.GetReconstructedBase(node.GetId(), site_idx);
      // Reconstructed base should always be unambiguous
      TestAssert(base == 'A' || base == 'C' || base == 'G' || base == 'T');
    }
  }
}

[[maybe_unused]] static const auto test_ambiguous_leaf_bases_registered =
    add_test({test_ambiguous_leaf_bases, "Sankoff: Ambiguous leaf bases", {"sankoff"}});

// Test that reconstruction must be called after scoring
static void test_reconstruction_order() {
  auto dag_storage = make_sample_dag();
  auto dag = dag_storage.View();

  SankoffScorer<decltype(dag)> scorer(dag);
  auto root = dag.GetRoot();

  // Should throw if calling ReconstructAncestralSequences before ComputeScoreBelow
  TestThrow(scorer.ReconstructAncestralSequences(root));

  // Should throw if calling GetReconstructedGenome before reconstruction
  scorer.ComputeScoreBelow(root);
  TestThrow(scorer.GetReconstructedGenome(root));

  // After proper sequence, should work
  scorer.ReconstructAncestralSequences(root);
  auto genome = scorer.GetReconstructedGenome(root);  // Should not throw
}

[[maybe_unused]] static const auto test_reconstruction_order_registered =
    add_test({test_reconstruction_order, "Sankoff: Reconstruction order", {"sankoff"}});

// Test with topology-only DAG (no mutations)
static void test_no_mutations() {
  auto dag_storage = make_sample_dag_topology();
  auto dag = dag_storage.View();

  SankoffScorer<decltype(dag)> scorer(dag);

  // No variable sites when there are no mutations
  TestAssert(scorer.GetNumVariableSites() == 0);

  // Scoring should return 0
  auto root = dag.GetRoot();
  double score = scorer.ComputeScoreBelow(root);
  TestAssert(score == 0.0);
}

[[maybe_unused]] static const auto test_no_mutations_registered =
    add_test({test_no_mutations, "Sankoff: No mutations", {"sankoff"}});

// Test comparing Sankoff uniform score with manual mutation count
static void test_score_at_most_mutation_count() {
  auto dag_storage = make_unambiguous_sample_dag();
  auto dag = dag_storage.View();

  SankoffScorer<decltype(dag)> scorer(dag);

  auto root = dag.GetRoot();
  double sankoff_score = scorer.ComputeScoreBelow(root);

  // Count mutations on all edges manually
  size_t mutation_count = 0;
  for (auto edge : dag.GetEdges()) {
    mutation_count += edge.GetEdgeMutations().size();
  }

  // Sankoff finds the optimal ancestral reconstruction, so its score should be
  // at most the mutation count from the given (possibly suboptimal) assignment.
  // It may be less if Sankoff finds a better reconstruction.
  TestAssert(sankoff_score <= static_cast<double>(mutation_count));
  TestAssert(sankoff_score > 0.0);  // Should have some cost
}

[[maybe_unused]] static const auto test_score_at_most_mutation_count_registered =
    add_test({test_score_at_most_mutation_count,
              "Sankoff: Score at most mutation count",
              {"sankoff"}});

// ============================================================================
// Validation tests against Python historydag ground truth
// ============================================================================

// Helper: build a 4-site, 3-leaf tree for validation tests.
//
// Topology:
//     UA(0) -> root(1) -> { n1(2) -> {L1(3), L2(4)}, L3(5) }
//
// Reference sequence: "AAAA"
// Leaf sequences: L1="CTAA", L2="CACC", L3="AGTC"
[[maybe_unused]] static auto make_four_site_validation_dag() {
  MADAGStorage<> storage = MADAGStorage<>::EmptyDefault();
  auto dag = storage.View();

  dag.SetReferenceSequence("AAAA");
  dag.InitializeNodes(6);

  size_t edge_id = 0;
  dag.AddEdge({edge_id++}, {0}, {1}, {0});  // UA(0) -> root(1)
  dag.AddEdge({edge_id++}, {1}, {2}, {0});  // root(1) -> n1(2)
  dag.AddEdge({edge_id++}, {1}, {5}, {1});  // root(1) -> L3(5)
  dag.AddEdge({edge_id++}, {2}, {3}, {0});  // n1(2) -> L1(3)
  dag.AddEdge({edge_id++}, {2}, {4}, {1});  // n1(2) -> L2(4)

  dag.BuildConnections();

  NodeSeqMap seq_map;
  seq_map[{0}] = {"AAAA"};  // UA
  seq_map[{1}] = {"AAAA"};  // root
  seq_map[{2}] = {"AAAA"};  // n1
  seq_map[{3}] = {"CTAA"};  // L1
  seq_map[{4}] = {"CACC"};  // L2
  seq_map[{5}] = {"AGTC"};  // L3

  dag.SetCompactGenomesFromNodeSequenceMap(seq_map);
  dag.RecomputeEdgeMutations();

  return storage;
}

// Test with 4 different per-site cost matrices.
// Ground truth verified with Python historydag package.
// Expected per-site scores: 2.0, 2.0, 6.0, 2.0 => total 12.0
static void test_sankoff_four_sites() {
  auto dag_storage = make_four_site_validation_dag();
  auto dag = dag_storage.View();

  // Site 1 (pos 1): transition/transversion — transitions cost 1, transversions 2
  SiteCostMatrix site1;
  site1[0] = {0.0, 2.0, 1.0, 2.0};  // from A
  site1[1] = {2.0, 0.0, 2.0, 1.0};  // from C
  site1[2] = {1.0, 2.0, 0.0, 2.0};  // from G
  site1[3] = {2.0, 1.0, 2.0, 0.0};  // from T

  // Site 2 (pos 2): uniform
  SiteCostMatrix site2;
  site2[0] = {0.0, 1.0, 1.0, 1.0};
  site2[1] = {1.0, 0.0, 1.0, 1.0};
  site2[2] = {1.0, 1.0, 0.0, 1.0};
  site2[3] = {1.0, 1.0, 1.0, 0.0};

  // Site 3 (pos 3): A<->G very cheap (0.5), everything else 3
  SiteCostMatrix site3;
  site3[0] = {0.0, 3.0, 0.5, 3.0};
  site3[1] = {3.0, 0.0, 3.0, 3.0};
  site3[2] = {0.5, 3.0, 0.0, 3.0};
  site3[3] = {3.0, 3.0, 3.0, 0.0};

  // Site 4 (pos 4): asymmetric — A->C=1 but C->A=2
  SiteCostMatrix site4;
  site4[0] = {0.0, 1.0, 1.0, 1.0};
  site4[1] = {2.0, 0.0, 1.0, 1.0};  // C->A = 2 (asymmetric)
  site4[2] = {1.0, 1.0, 0.0, 1.0};
  site4[3] = {1.0, 1.0, 1.0, 0.0};

  std::unordered_map<size_t, SiteCostMatrix> site_matrices;
  site_matrices[1] = site1;
  site_matrices[2] = site2;
  site_matrices[3] = site3;
  site_matrices[4] = site4;

  SankoffScorer<decltype(dag)> scorer(dag, site_matrices);

  auto root = dag.GetRoot();
  double total_score = scorer.ComputeScoreBelow(root);

  // Verify total score matches Python ground truth
  TestAssert(std::abs(total_score - 12.0) < 1e-9);

  // Verify per-site scores from DP table
  const auto& dp_table = scorer.GetDPTable();
  size_t root_idx = root.GetId().value;

  double expected_site_scores[] = {2.0, 2.0, 6.0, 2.0};
  double sum_of_sites = 0.0;

  for (size_t site_idx = 0; site_idx < 4; ++site_idx) {
    const SiteCosts& costs = dp_table.dp_costs[root_idx][site_idx];
    double min_cost = *std::min_element(costs.begin(), costs.end());
    sum_of_sites += min_cost;
    TestAssert(std::abs(min_cost - expected_site_scores[site_idx]) < 1e-9);
  }

  TestAssert(std::abs(sum_of_sites - 12.0) < 1e-9);

  // Verify leaf bases are preserved after reconstruction
  scorer.ReconstructAncestralSequences(root);

  std::unordered_map<size_t, std::string> expected_leaf_seqs = {
      {3, "CTAA"}, {4, "CACC"}, {5, "AGTC"}};

  for (const auto& [node_id, expected_seq] : expected_leaf_seqs) {
    for (size_t site_idx = 0; site_idx < 4; ++site_idx) {
      char recon = scorer.GetReconstructedBase(NodeId{node_id}, site_idx);
      TestAssert(recon == expected_seq[site_idx]);
    }
  }
}

[[maybe_unused]] static const auto test_sankoff_four_sites_registered =
    add_test({test_sankoff_four_sites,
              "Sankoff: Four sites with different matrices",
              {"sankoff"}});

// Test that Sankoff with uniform costs matches hand-calculated Fitch optimal.
// Same 4-site, 3-leaf tree.
//
// With uniform costs (0 diagonal, 1 off-diagonal), Sankoff should find the
// Fitch-optimal score. Hand calculation:
//   Site 1 (C,C,A): n1=C(0), root chooses A or C -> 1 change minimum => 1
//   Site 2 (T,A,G): all different, n1 needs 1 change, root needs 1 more => 2
//   Site 3 (A,C,T): all different, n1 needs 1 change, root needs 1 more => 2
//   Site 4 (A,C,C): n1 has {A,C} cost 1; or C cost 0+0=0 from L2,
//                    actually n1=C(1), root=C(0+1=1) => 1
//   Total Fitch-optimal = 1 + 2 + 2 + 1 = 6
static void test_sankoff_uniform_matches_fitch() {
  auto dag_storage = make_four_site_validation_dag();
  auto dag = dag_storage.View();

  SankoffScorer<decltype(dag)> scorer(dag);

  auto root = dag.GetRoot();
  double sankoff_score = scorer.ComputeScoreBelow(root);

  // Sankoff with uniform costs should give Fitch-optimal score
  TestAssert(std::abs(sankoff_score - 6.0) < 1e-9);

  // Also verify: Sankoff score <= ParsimonyScore from edge annotations.
  // The edge mutations were computed from internal nodes set to reference "AAAA",
  // so ParsimonyScore counts all mutations on edges (likely suboptimal).
  MADAG const_dag = dag;
  SubtreeWeight<ParsimonyScore, MADAG> weight{const_dag};
  ParsimonyScore::Weight parsimony =
      weight.ComputeWeightBelow(const_dag.GetRoot(), {});
  TestAssert(sankoff_score <= static_cast<double>(parsimony));
}

[[maybe_unused]] static const auto test_sankoff_uniform_matches_fitch_registered =
    add_test({test_sankoff_uniform_matches_fitch,
              "Sankoff: Uniform costs match Fitch optimal",
              {"sankoff"}});

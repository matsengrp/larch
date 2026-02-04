#include "larch/subtree/sankoff.hpp"

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

  // With uniform costs, the score should equal the parsimony score
  // (number of mutations on the optimal tree)
  // The sample DAG is actually a tree (single path from root to each leaf)
  // We can verify the score is non-negative and reasonable
  TestAssert(score >= 0.0);
  TestAssert(score < 1000.0);  // Sanity check

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

[[maybe_unused]] static const auto test_leaf_base_assignment_registered = add_test(
    {test_leaf_base_assignment, "Sankoff: Leaf base assignment", {"sankoff"}});

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
        ts_tv_matrix[i][j] = 1.0;  // Transition
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
              "Sankoff: Transition/transversion costs", {"sankoff"}});

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

[[maybe_unused]] static const auto test_site_specific_costs_registered = add_test(
    {test_site_specific_costs, "Sankoff: Site-specific costs", {"sankoff"}});

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

[[maybe_unused]] static const auto test_ambiguous_leaf_bases_registered = add_test(
    {test_ambiguous_leaf_bases, "Sankoff: Ambiguous leaf bases", {"sankoff"}});

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

[[maybe_unused]] static const auto test_reconstruction_order_registered = add_test(
    {test_reconstruction_order, "Sankoff: Reconstruction order", {"sankoff"}});

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
              "Sankoff: Score at most mutation count", {"sankoff"}});

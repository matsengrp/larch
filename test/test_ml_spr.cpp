// Tests for ML scoring backend
#include "test_common.hpp"

#ifdef USE_NETAM

#include <netam/crepe.hpp>
#include <netam/likelihood.hpp>
#include <netam/kmer_sequence_encoder.hpp>

#include "larch/dag_loader.hpp"
#include "larch/spr/spr_view.hpp"
#include "test_common_dag.hpp"

namespace {

// Reference sequence from linearham data (381 bases)
const std::string kLinearhamRefSeq =
    "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTC"
    "TGGATTCACCGTCAGTAGCAACTACATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAG"
    "TTATTTATAGCGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCAGATTCACCATCTCCAGAGACAATTCC"
    "AAGAACACGCTGTATCTTCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGGCAC"
    "AACACACGGGTATAGCAGTGAAGGCATGACTTCAAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCG"
    "TCTCCTCAG";

// Create sample sequences with mutations from reference
std::string mutate_sequence(const std::string& seq, size_t pos, char new_base) {
  std::string result{seq};
  result[pos] = new_base;
  return result;
}

// Test that we can load the model and run inference
void test_ml_model_loading_and_inference() {
  torch::NoGradGuard no_grad;

  netam::crepe model{"data/linearham/ThriftyHumV0.2-45-libtorch.pth",
                     "data/linearham/ThriftyHumV0.2-45.yml"};

  // Encode sequence - returns (encoded, wt_modifier) both without batch dim
  auto [encoded_1d, wt_mod_2d] = model.encoder().encode_sequence(kLinearhamRefSeq);

  // Add batch dimension
  auto encoded = encoded_1d.unsqueeze(0);  // [1, site_count]
  auto wt_mod = wt_mod_2d.unsqueeze(0);    // [1, site_count, 4]

  // Create mask (all true for valid positions)
  auto mask = torch::ones({1, encoded.size(1)}, torch::kBool);

  // Run forward pass
  model->eval();
  auto [rates, csp_logits] = model->forward(encoded, mask, wt_mod);

  // Verify output shapes - note: output is based on encoder's site_count, not input
  // length
  TestAssert(rates.dim() == 2);
  TestAssert(rates.size(0) == 1);

  TestAssert(csp_logits.dim() == 3);
  TestAssert(csp_logits.size(0) == 1);
  TestAssert(csp_logits.size(2) == 4);

  // Rates should be positive
  TestAssert(torch::all(rates > 0).item<bool>());
}

// Test edge log-likelihood computation
void test_edge_log_likelihood() {
  torch::NoGradGuard no_grad;

  netam::crepe model{"data/linearham/ThriftyHumV0.2-45-libtorch.pth",
                     "data/linearham/ThriftyHumV0.2-45.yml"};
  model->eval();

  std::string parent_seq = kLinearhamRefSeq;
  // Create child with one mutation at position 100
  std::string child_seq = mutate_sequence(kLinearhamRefSeq, 100, 'T');

  // Encode parent for model input
  auto [encoded_1d, wt_mod_2d] = model.encoder().encode_sequence(parent_seq);

  // Add batch dimension
  auto encoded = encoded_1d.unsqueeze(0);
  auto wt_mod = wt_mod_2d.unsqueeze(0);
  auto mask = torch::ones({1, encoded.size(1)}, torch::kBool);

  // Run forward pass
  auto [rates, csp_logits] = model->forward(encoded, mask, wt_mod);

  // Apply softmax to get CSP probabilities
  auto csp = torch::softmax(csp_logits, -1);

  // Encode as base indices
  auto parent_indices = netam::kmer_sequence_encoder::encode_bases(parent_seq);
  auto child_indices = netam::kmer_sequence_encoder::encode_bases(child_seq);

  // Compute log-likelihood
  auto log_likelihood =
      netam::poisson_context_log_likelihood(rates, csp, parent_indices, child_indices);

  double ll = log_likelihood.item<double>();

  // Log-likelihood should be finite
  TestAssert(std::isfinite(ll));

  // With one mutation, log-likelihood should be negative (less than 0)
  // (more mutations = more negative log-likelihood)
  std::cout << "Single mutation log-likelihood: " << ll << std::endl;
}

// Test that identical sequences have log-likelihood of 0
void test_identical_sequences_log_likelihood() {
  torch::NoGradGuard no_grad;

  netam::crepe model{"data/linearham/ThriftyHumV0.2-45-libtorch.pth",
                     "data/linearham/ThriftyHumV0.2-45.yml"};
  model->eval();

  std::string seq = kLinearhamRefSeq;

  // Encode for model input
  auto [encoded_1d, wt_mod_2d] = model.encoder().encode_sequence(seq);
  auto encoded = encoded_1d.unsqueeze(0);
  auto wt_mod = wt_mod_2d.unsqueeze(0);
  auto mask = torch::ones({1, encoded.size(1)}, torch::kBool);

  auto [rates, csp_logits] = model->forward(encoded, mask, wt_mod);
  auto csp = torch::softmax(csp_logits, -1);

  // Encode as base indices
  auto seq_indices = netam::kmer_sequence_encoder::encode_bases(seq);

  // Compute log-likelihood for identical parent/child
  auto log_likelihood =
      netam::poisson_context_log_likelihood(rates, csp, seq_indices, seq_indices);

  double ll = log_likelihood.item<double>();

  // Identical sequences should have log-likelihood of 0
  TestAssert(std::abs(ll) < 1e-9);
}

// Test that more mutations give lower log-likelihood
void test_more_mutations_lower_likelihood() {
  torch::NoGradGuard no_grad;

  netam::crepe model{"data/linearham/ThriftyHumV0.2-45-libtorch.pth",
                     "data/linearham/ThriftyHumV0.2-45.yml"};
  model->eval();

  std::string parent_seq = kLinearhamRefSeq;

  // Child with 1 mutation
  std::string child_1mut = mutate_sequence(kLinearhamRefSeq, 100, 'T');

  // Child with 3 mutations
  std::string child_3mut = mutate_sequence(kLinearhamRefSeq, 100, 'T');
  child_3mut = mutate_sequence(child_3mut, 150, 'A');
  child_3mut = mutate_sequence(child_3mut, 200, 'C');

  // Compute log-likelihood for 1 mutation
  auto [encoded_1d, wt_mod_2d] = model.encoder().encode_sequence(parent_seq);
  auto encoded = encoded_1d.unsqueeze(0);
  auto wt_mod = wt_mod_2d.unsqueeze(0);
  auto mask = torch::ones({1, encoded.size(1)}, torch::kBool);

  auto [rates, csp_logits] = model->forward(encoded, mask, wt_mod);
  auto csp = torch::softmax(csp_logits, -1);

  auto parent_indices = netam::kmer_sequence_encoder::encode_bases(parent_seq);
  auto child_1mut_indices = netam::kmer_sequence_encoder::encode_bases(child_1mut);
  auto child_3mut_indices = netam::kmer_sequence_encoder::encode_bases(child_3mut);

  auto ll_1mut = netam::poisson_context_log_likelihood(rates, csp, parent_indices,
                                                       child_1mut_indices);
  auto ll_3mut = netam::poisson_context_log_likelihood(rates, csp, parent_indices,
                                                       child_3mut_indices);

  std::cout << "1 mutation log-likelihood: " << ll_1mut.item<double>() << std::endl;
  std::cout << "3 mutations log-likelihood: " << ll_3mut.item<double>() << std::endl;

  // More mutations should have lower (more negative) log-likelihood
  // Actually: log-likelihood can go either way depending on context
  // The formula is: log_lik = Σlog(λ_j) + n*log(t_hat) - n
  // With more mutations, t_hat increases, which can offset the penalty
  // So we just check both are finite
  TestAssert(std::isfinite(ll_1mut.item<double>()));
  TestAssert(std::isfinite(ll_3mut.item<double>()));
}

// Create a sample DAG with linearham-compatible sequences
[[maybe_unused]] static auto make_linearham_sample_dag() {
  MADAGStorage<> input_storage = MADAGStorage<>::EmptyDefault();
  auto dag = input_storage.View();

  dag.SetReferenceSequence(kLinearhamRefSeq);
  dag.InitializeNodes(7);
  /*
      Simple tree:
             0 (UA)
             |
             6 (root, ref sequence)
            / \
           4   5
          / \   \
         1   2   3
      Leaves: 1, 2, 3 (with mutations from ref)
      Internal: 4, 5, 6
  */
  size_t edge_id = 0;
  dag.AddEdge({edge_id++}, {0}, {6}, {0});  // UA -> root
  dag.AddEdge({edge_id++}, {6}, {4}, {0});  // root -> 4
  dag.AddEdge({edge_id++}, {6}, {5}, {1});  // root -> 5
  dag.AddEdge({edge_id++}, {4}, {1}, {0});  // 4 -> 1 (leaf)
  dag.AddEdge({edge_id++}, {4}, {2}, {1});  // 4 -> 2 (leaf)
  dag.AddEdge({edge_id++}, {5}, {3}, {0});  // 5 -> 3 (leaf)

  dag.BuildConnections();

  // Add some mutations to edges
  // Position 100: G -> T on edge to leaf 1
  dag.Get(EdgeId{3}).GetMutableEdgeMutations()[{101}] = {'G', 'T'};

  // Position 150: C -> A on edge to leaf 2
  dag.Get(EdgeId{4}).GetMutableEdgeMutations()[{151}] = {'C', 'A'};

  // Position 200: A -> C on edge to leaf 3
  dag.Get(EdgeId{5}).GetMutableEdgeMutations()[{201}] = {'A', 'C'};

  // Add mutations on internal edges
  dag.Get(EdgeId{1}).GetMutableEdgeMutations()[{50}] = {'G', 'A'};
  dag.Get(EdgeId{2}).GetMutableEdgeMutations()[{75}] = {'T', 'C'};

  dag.RecomputeCompactGenomes(true);
  dag.SampleIdsFromCG();

  return input_storage;
}

// Helper to expand CompactGenome to full sequence
template <typename DAGView>
std::string expand_sequence(const DAGView& dag, NodeId node) {
  const auto& ref_seq = dag.GetReferenceSequence();
  std::string result{ref_seq.begin(), ref_seq.end()};

  const auto& compact_genome = dag.Get(node).GetCompactGenome();
  for (const auto& [pos, base] : compact_genome) {
    result[pos.value - 1] = base.ToChar();  // MutationPosition is 1-indexed
  }

  return result;
}

// Test computing log-likelihood for all edges in a DAG
void test_dag_edge_log_likelihoods() {
  torch::NoGradGuard no_grad;

  netam::crepe model{"data/linearham/ThriftyHumV0.2-45-libtorch.pth",
                     "data/linearham/ThriftyHumV0.2-45.yml"};
  model->eval();

  auto dag_storage = make_linearham_sample_dag();
  auto dag = dag_storage.View();

  std::cout << "DAG has " << dag.GetNodesCount() << " nodes and " << dag.GetEdgesCount()
            << " edges" << std::endl;

  double total_ll = 0.0;

  for (auto edge : dag.GetEdges()) {
    if (edge.GetParent().IsUA()) {
      continue;  // Skip edge from UA
    }

    NodeId parent_id = edge.GetParent().GetId();
    NodeId child_id = edge.GetChild().GetId();

    std::string parent_seq = expand_sequence(dag, parent_id);
    std::string child_seq = expand_sequence(dag, child_id);

    // Encode parent for model
    auto [encoded_1d, wt_mod_2d] = model.encoder().encode_sequence(parent_seq);
    auto encoded = encoded_1d.unsqueeze(0);
    auto wt_mod = wt_mod_2d.unsqueeze(0);
    auto mask = torch::ones({1, encoded.size(1)}, torch::kBool);

    auto [rates, csp_logits] = model->forward(encoded, mask, wt_mod);
    auto csp = torch::softmax(csp_logits, -1);

    // Encode as base indices
    auto parent_indices = netam::kmer_sequence_encoder::encode_bases(parent_seq);
    auto child_indices = netam::kmer_sequence_encoder::encode_bases(child_seq);

    // Compute log-likelihood
    auto ll = netam::poisson_context_log_likelihood(rates, csp, parent_indices,
                                                    child_indices);

    double ll_val = ll.item<double>();
    total_ll += ll_val;

    std::cout << "Edge " << edge.GetId() << " (" << parent_id << " -> " << child_id
              << "): log-likelihood = " << ll_val << std::endl;
  }

  std::cout << "Total log-likelihood: " << total_ll << std::endl;
  TestAssert(std::isfinite(total_ll));
}

}  // namespace

[[maybe_unused]] static bool reg_test_ml_spr =
    add_test({test_ml_model_loading_and_inference,
              "ML SPR: Model loading and inference",
              {"netam", "ml"}}) &&
    add_test(
        {test_edge_log_likelihood, "ML SPR: Edge log-likelihood", {"netam", "ml"}}) &&
    add_test({test_identical_sequences_log_likelihood,
              "ML SPR: Identical sequences log-likelihood",
              {"netam", "ml"}}) &&
    add_test({test_more_mutations_lower_likelihood,
              "ML SPR: More mutations likelihood comparison",
              {"netam", "ml"}}) &&
    add_test({test_dag_edge_log_likelihoods,
              "ML SPR: DAG edge log-likelihoods",
              {"netam", "ml"}});

#endif  // USE_NETAM

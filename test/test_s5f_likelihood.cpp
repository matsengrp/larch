// Tests for S5F model log-likelihood computation
// Compares C++ implementation against Python ground truth values

#include "test_common.hpp"

#ifdef USE_NETAM

#include <netam/crepe.hpp>
#include <netam/likelihood.hpp>
#include <netam/kmer_sequence_encoder.hpp>

#include <cmath>

namespace {

// Reference sequence (379 bases)
const std::string kReferenceSeq =
    "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTC"
    "TGGATTCACCGTCAGTAGCAACTACATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAG"
    "TTATTTATAGCGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCAGATTCACCATCTCCAGAGACAATTCC"
    "AAGAACACGCTGTATCTTCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGGCAC"
    "AACACACGGGTATAGCAGTGAAGGCATGACTTCAAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCG"
    "TCTCCTCAG";

// Ground truth log-likelihoods computed with S5F model from Python
// These values should match the C++ implementation
struct S5FGroundTruth {
  const char* name;
  std::size_t mutation_pos;
  char old_base;
  char new_base;
  double log_likelihood;
};

// Ground truth - note: position indices are 0-based
// The sequence at positions 10, 50, 100, 150, 200 are: T, C, T, A, C
const S5FGroundTruth kS5FGroundTruth[] = {
    {"edge_1_T_to_A_pos10", 10, 'T', 'A', -8.55931615829468},
    {"edge_2_C_to_T_pos50", 50, 'C', 'T', -7.91756159067154},
    {"edge_3_T_to_G_pos100", 100, 'T', 'G', -9.04228067398071},
    {"edge_4_A_to_C_pos150", 150, 'A', 'C', -7.1298691034317},
    {"edge_5_C_to_A_pos200", 200, 'C', 'A', -6.76613847911358},
};

// Create a sequence with a single mutation
std::string mutate_sequence(const std::string& seq, std::size_t pos, char new_base) {
  std::string result{seq};
  result[pos] = new_base;
  return result;
}

// Compute log-likelihood for a parent-child pair using the Poisson context model
double compute_log_likelihood(netam::crepe& model, const std::string& parent_seq,
                              const std::string& child_seq) {
  torch::NoGradGuard no_grad;

  // Encode parent sequence
  auto [encoded_1d, wt_mod_2d] = model.encoder().encode_sequence(parent_seq);

  // Add batch dimension
  auto encoded = encoded_1d.unsqueeze(0);
  auto wt_mod = wt_mod_2d.unsqueeze(0);
  auto mask = torch::ones({1, encoded.size(1)}, torch::kBool);

  // Run forward pass
  model->eval();
  auto [rates, csp_logits] = model->forward(encoded, mask, wt_mod);

  // Apply softmax to get CSP probabilities
  auto csp = torch::softmax(csp_logits, -1);

  // Encode as base indices
  auto parent_indices = netam::kmer_sequence_encoder::encode_bases(parent_seq);
  auto child_indices = netam::kmer_sequence_encoder::encode_bases(child_seq);

  // Compute log-likelihood
  auto log_likelihood =
      netam::poisson_context_log_likelihood(rates, csp, parent_indices, child_indices);

  return log_likelihood.item<double>();
}

// Test that S5F model can be loaded
void test_s5f_model_loading() {
  torch::NoGradGuard no_grad;

  netam::crepe model{"data/linearham/s5f-libtorch.pth", "data/linearham/s5f.yml"};

  // Encode sequence
  auto [encoded_1d, wt_mod_2d] = model.encoder().encode_sequence(kReferenceSeq);

  // Add batch dimension
  auto encoded = encoded_1d.unsqueeze(0);
  auto wt_mod = wt_mod_2d.unsqueeze(0);
  auto mask = torch::ones({1, encoded.size(1)}, torch::kBool);

  // Run forward pass
  model->eval();
  auto [rates, csp_logits] = model->forward(encoded, mask, wt_mod);

  // Verify output shapes
  TestAssert(rates.dim() == 2);
  TestAssert(rates.size(0) == 1);

  TestAssert(csp_logits.dim() == 3);
  TestAssert(csp_logits.size(0) == 1);
  TestAssert(csp_logits.size(2) == 4);

  // Rates should be positive
  TestAssert(torch::all(rates > 0).item<bool>());

  std::cout << "S5F model loaded successfully" << std::endl;
  std::cout << "  Output rates shape: [" << rates.size(0) << ", " << rates.size(1)
            << "]" << std::endl;
  std::cout << "  CSP logits shape: [" << csp_logits.size(0) << ", "
            << csp_logits.size(1) << ", " << csp_logits.size(2) << "]" << std::endl;
}

// Test identical sequences have zero log-likelihood
void test_s5f_identical_sequences() {
  torch::NoGradGuard no_grad;

  netam::crepe model{"data/linearham/s5f-libtorch.pth", "data/linearham/s5f.yml"};

  double ll = compute_log_likelihood(model, kReferenceSeq, kReferenceSeq);

  std::cout << "Identical sequences log-likelihood: " << ll << std::endl;
  TestAssert(std::abs(ll) < 1e-9);
}

// Test single mutations against Python ground truth
void test_s5f_single_mutations() {
  torch::NoGradGuard no_grad;

  netam::crepe model{"data/linearham/s5f-libtorch.pth", "data/linearham/s5f.yml"};

  const double tolerance = 1e-5;  // Allow small floating point differences

  std::cout << "Testing single mutations against Python ground truth:" << std::endl;

  for (const auto& truth : kS5FGroundTruth) {
    // Verify the original base is what we expect
    TestAssert(kReferenceSeq[truth.mutation_pos] == truth.old_base);

    std::string child_seq =
        mutate_sequence(kReferenceSeq, truth.mutation_pos, truth.new_base);

    double ll = compute_log_likelihood(model, kReferenceSeq, child_seq);

    double diff = std::abs(ll - truth.log_likelihood);

    std::cout << "  " << truth.name << ": C++=" << ll
              << " Python=" << truth.log_likelihood << " diff=" << diff << std::endl;

    if (diff > tolerance) {
      std::cerr << "MISMATCH: " << truth.name << " - C++ result " << ll
                << " differs from Python ground truth " << truth.log_likelihood
                << " by " << diff << std::endl;
    }
    TestAssert(diff < tolerance);
  }
}

// Test multiple mutations
void test_s5f_multiple_mutations() {
  torch::NoGradGuard no_grad;

  netam::crepe model{"data/linearham/s5f-libtorch.pth", "data/linearham/s5f.yml"};

  // Multiple mutations at positions 10, 50, 100
  // Positions: 10=T, 50=C, 100=T
  std::string multi_mut = mutate_sequence(kReferenceSeq, 10, 'G');  // T->G
  multi_mut = mutate_sequence(multi_mut, 50, 'A');                  // C->A
  multi_mut = mutate_sequence(multi_mut, 100, 'C');                 // T->C

  double ll = compute_log_likelihood(model, kReferenceSeq, multi_mut);

  // Python ground truth: -22.4804043769836
  const double expected = -22.4804043769836;
  const double tolerance = 1e-5;

  double diff = std::abs(ll - expected);
  std::cout << "Multiple mutations: C++=" << ll << " Python=" << expected
            << " diff=" << diff << std::endl;

  TestAssert(diff < tolerance);
}

}  // namespace

[[maybe_unused]] static bool reg_test_s5f =
    add_test({test_s5f_model_loading, "S5F: Model loading", {"netam", "s5f"}}) &&
    add_test({test_s5f_identical_sequences,
              "S5F: Identical sequences log-likelihood",
              {"netam", "s5f"}}) &&
    add_test({test_s5f_single_mutations,
              "S5F: Single mutations vs Python ground truth",
              {"netam", "s5f"}}) &&
    add_test({test_s5f_multiple_mutations,
              "S5F: Multiple mutations vs Python ground truth",
              {"netam", "s5f"}});

#endif  // USE_NETAM

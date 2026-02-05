// Tests for indep_rscnn_model
#include "test_common.hpp"

#ifdef USE_NETAM
#include <netam/indep_rscnn_model.hpp>

#include <torch/torch.h>

using namespace netam;

namespace {

// Helper to create a minimal YAML config for testing
YAML::Node make_config(std::size_t kmer_length = 3, std::size_t embedding_dim = 7,
                       std::size_t filter_count = 16, std::size_t kernel_size = 9,
                       double dropout_prob = 0.2) {
  YAML::Node yaml;
  yaml["kmer_length"] = kmer_length;
  yaml["embedding_dim"] = embedding_dim;
  yaml["filter_count"] = filter_count;
  yaml["kernel_size"] = kernel_size;
  yaml["dropout_prob"] = dropout_prob;
  return yaml;
}

// ============================================================================
// indep_rscnn_params tests
// ============================================================================

void test_params_construction() {
  torch::NoGradGuard no_grad;
  auto yaml = make_config(3, 7, 16, 9, 0.2);
  indep_rscnn_params params{65, yaml};

  TestAssert(params.kmer_count() == 65);
  TestAssert(params.kmer_length() == 3);
  TestAssert(params.embedding_dim() == 7);
  TestAssert(params.filter_count() == 16);
  TestAssert(params.kernel_size() == 9);
  TestAssert(std::abs(params.dropout_prob() - 0.2) < 1e-9);
}

void test_params_different_values() {
  torch::NoGradGuard no_grad;
  auto yaml = make_config(5, 32, 64, 15, 0.5);
  indep_rscnn_params params{1025, yaml};

  TestAssert(params.kmer_count() == 1025);
  TestAssert(params.kmer_length() == 5);
  TestAssert(params.embedding_dim() == 32);
  TestAssert(params.filter_count() == 64);
  TestAssert(params.kernel_size() == 15);
  TestAssert(std::abs(params.dropout_prob() - 0.5) < 1e-9);
}

// ============================================================================
// indep_rscnn_model tests
// ============================================================================

void test_model_output_shapes() {
  torch::NoGradGuard no_grad;
  auto yaml = make_config();
  indep_rscnn_params params{65, yaml};
  indep_rscnn_model model{params};
  model.eval();

  const int64_t batch_size = 2;
  const int64_t seq_length = 100;

  auto encoded = torch::randint(0, 65, {batch_size, seq_length}, torch::kInt32);
  auto mask = torch::ones({batch_size, seq_length}, torch::kBool);
  auto wt_modifier = torch::zeros({batch_size, seq_length, 4});

  auto [rates, csp_logits] = model.forward(encoded, mask, wt_modifier);

  // rates should be [B, L]
  TestAssert(rates.dim() == 2);
  TestAssert(rates.size(0) == batch_size);
  TestAssert(rates.size(1) == seq_length);

  // csp_logits should be [B, L, 4]
  TestAssert(csp_logits.dim() == 3);
  TestAssert(csp_logits.size(0) == batch_size);
  TestAssert(csp_logits.size(1) == seq_length);
  TestAssert(csp_logits.size(2) == 4);
}

void test_model_batch_size_1() {
  torch::NoGradGuard no_grad;
  auto yaml = make_config();
  indep_rscnn_params params{65, yaml};
  indep_rscnn_model model{params};
  model.eval();

  const int64_t batch_size = 1;
  const int64_t seq_length = 50;

  auto encoded = torch::randint(0, 65, {batch_size, seq_length}, torch::kInt32);
  auto mask = torch::ones({batch_size, seq_length}, torch::kBool);
  auto wt_modifier = torch::zeros({batch_size, seq_length, 4});

  auto [rates, csp_logits] = model.forward(encoded, mask, wt_modifier);

  TestAssert(rates.size(0) == 1);
  TestAssert(rates.size(1) == 50);
  TestAssert(csp_logits.size(0) == 1);
  TestAssert(csp_logits.size(1) == 50);
}

void test_model_rates_positive() {
  torch::NoGradGuard no_grad;
  auto yaml = make_config();
  indep_rscnn_params params{65, yaml};
  indep_rscnn_model model{params};
  model.eval();

  auto encoded = torch::randint(0, 65, {1, 100}, torch::kInt32);
  auto mask = torch::ones({1, 100}, torch::kBool);
  auto wt_modifier = torch::zeros({1, 100, 4});

  auto [rates, csp_logits] = model.forward(encoded, mask, wt_modifier);

  // Rates are exp(log_rates), so should all be positive
  TestAssert(torch::all(rates > 0).item<bool>());
}

void test_model_mask_zeros_output() {
  torch::NoGradGuard no_grad;
  auto yaml = make_config();
  indep_rscnn_params params{65, yaml};
  indep_rscnn_model model{params};
  model.eval();

  auto encoded = torch::randint(0, 65, {1, 10}, torch::kInt32);
  auto wt_modifier = torch::zeros({1, 10, 4});

  // Mask with some positions set to false
  auto mask = torch::ones({1, 10}, torch::kBool);
  mask[0][5] = false;
  mask[0][6] = false;

  auto [rates, csp_logits] = model.forward(encoded, mask, wt_modifier);

  // Where mask is false (0), rates should be exp(0) = 1
  TestAssert(std::abs(rates[0][5].item<float>() - 1.0f) < 1e-5f);
  TestAssert(std::abs(rates[0][6].item<float>() - 1.0f) < 1e-5f);

  // csp_logits should be zero where masked
  TestAssert(torch::all(csp_logits[0][5] == torch::tensor(0)).item<bool>());
  TestAssert(torch::all(csp_logits[0][6] == torch::tensor(0)).item<bool>());
}

void test_model_wt_modifier_applied() {
  torch::NoGradGuard no_grad;
  auto yaml = make_config();
  indep_rscnn_params params{65, yaml};
  indep_rscnn_model model{params};
  model.eval();

  auto encoded = torch::randint(0, 65, {1, 10}, torch::kInt32);
  auto mask = torch::ones({1, 10}, torch::kBool);

  // Run without modifier
  auto wt_modifier_zero = torch::zeros({1, 10, 4});
  auto [rates1, csp_logits1] = model.forward(encoded, mask, wt_modifier_zero);

  // Run with modifier (large negative value at position 0, base 0)
  auto wt_modifier = torch::zeros({1, 10, 4});
  wt_modifier[0][0][0] = -1e9f;
  auto [rates2, csp_logits2] = model.forward(encoded, mask, wt_modifier);

  // Rates should be the same (wt_modifier only affects csp_logits)
  TestAssert(torch::allclose(rates1, rates2));

  // csp_logits at position 0, base 0 should be much lower
  TestAssert(csp_logits2[0][0][0].item<float>() < csp_logits1[0][0][0].item<float>());
  TestAssert(csp_logits2[0][0][0].item<float>() < -1e8f);
}

void test_model_adjust_rate_bias() {
  torch::NoGradGuard no_grad;
  auto yaml = make_config();
  indep_rscnn_params params{65, yaml};
  indep_rscnn_model model{params};
  model.eval();

  auto encoded = torch::randint(0, 65, {1, 10}, torch::kInt32);
  auto mask = torch::ones({1, 10}, torch::kBool);
  auto wt_modifier = torch::zeros({1, 10, 4});

  // Get rates before adjustment
  auto [rates_before, _1] = model.forward(encoded, mask, wt_modifier);

  // Adjust by log(2), which should double the rates
  model.adjust_rate_bias_by(std::log(2.0));

  // Get rates after adjustment
  auto [rates_after, _2] = model.forward(encoded, mask, wt_modifier);

  // rates_after should be approximately 2 * rates_before
  auto ratio = rates_after / rates_before;
  TestAssert(torch::allclose(ratio, torch::full_like(ratio, 2.0f),
                             /*rtol=*/1e-4, /*atol=*/1e-4));
}

void test_model_deterministic_in_eval_mode() {
  torch::NoGradGuard no_grad;
  auto yaml = make_config();
  indep_rscnn_params params{65, yaml};
  indep_rscnn_model model{params};
  model.eval();

  auto encoded = torch::randint(0, 65, {1, 50}, torch::kInt32);
  auto mask = torch::ones({1, 50}, torch::kBool);
  auto wt_modifier = torch::zeros({1, 50, 4});

  // Run twice, should get identical results in eval mode
  auto [rates1, csp_logits1] = model.forward(encoded, mask, wt_modifier);
  auto [rates2, csp_logits2] = model.forward(encoded, mask, wt_modifier);

  TestAssert(torch::equal(rates1, rates2));
  TestAssert(torch::equal(csp_logits1, csp_logits2));
}

void test_model_different_kernel_sizes() {
  torch::NoGradGuard no_grad;
  // Test with different kernel sizes to ensure padding works correctly
  for (std::size_t kernel_size : {3uz, 5uz, 7uz, 9uz, 11uz}) {
    auto yaml = make_config(3, 7, 16, kernel_size, 0.0);
    indep_rscnn_params params{65, yaml};
    indep_rscnn_model model{params};
    model.eval();

    auto encoded = torch::randint(0, 65, {1, 100}, torch::kInt32);
    auto mask = torch::ones({1, 100}, torch::kBool);
    auto wt_modifier = torch::zeros({1, 100, 4});

    auto [rates, csp_logits] = model.forward(encoded, mask, wt_modifier);

    // Output should maintain sequence length (same padding)
    TestAssert(rates.size(1) == 100);
    TestAssert(csp_logits.size(1) == 100);
  }
}

void test_model_short_sequence() {
  torch::NoGradGuard no_grad;
  auto yaml = make_config(3, 7, 16, 9, 0.0);
  indep_rscnn_params params{65, yaml};
  indep_rscnn_model model{params};
  model.eval();

  // Sequence shorter than kernel size
  const int64_t seq_length = 5;
  auto encoded = torch::randint(0, 65, {1, seq_length}, torch::kInt32);
  auto mask = torch::ones({1, seq_length}, torch::kBool);
  auto wt_modifier = torch::zeros({1, seq_length, 4});

  auto [rates, csp_logits] = model.forward(encoded, mask, wt_modifier);

  TestAssert(rates.size(1) == seq_length);
  TestAssert(csp_logits.size(1) == seq_length);
}

void test_model_csp_logits_sum_behavior() {
  torch::NoGradGuard no_grad;
  // After softmax, CSPs should sum to 1 for each position
  auto yaml = make_config();
  indep_rscnn_params params{65, yaml};
  indep_rscnn_model model{params};
  model.eval();

  auto encoded = torch::randint(0, 65, {1, 20}, torch::kInt32);
  auto mask = torch::ones({1, 20}, torch::kBool);
  auto wt_modifier = torch::zeros({1, 20, 4});

  auto [rates, csp_logits] = model.forward(encoded, mask, wt_modifier);

  // Apply softmax to get probabilities
  auto csps = torch::softmax(csp_logits, /*dim=*/-1);

  // Each position should sum to 1
  auto sums = csps.sum(/*dim=*/-1);
  TestAssert(torch::allclose(sums, torch::ones_like(sums), /*rtol=*/1e-5,
                             /*atol=*/1e-5));
}

void test_model_larger_batch() {
  torch::NoGradGuard no_grad;
  auto yaml = make_config();
  indep_rscnn_params params{65, yaml};
  indep_rscnn_model model{params};
  model.eval();

  const int64_t batch_size = 16;
  const int64_t seq_length = 100;

  auto encoded = torch::randint(0, 65, {batch_size, seq_length}, torch::kInt32);
  auto mask = torch::ones({batch_size, seq_length}, torch::kBool);
  auto wt_modifier = torch::zeros({batch_size, seq_length, 4});

  auto [rates, csp_logits] = model.forward(encoded, mask, wt_modifier);

  TestAssert(rates.size(0) == batch_size);
  TestAssert(csp_logits.size(0) == batch_size);
}

}  // namespace

[[maybe_unused]] static bool reg_test_model =
    add_test(
        {test_params_construction, "Netam Model: params construction", {"netam"}}) &&
    add_test({test_params_different_values,
              "Netam Model: params different values",
              {"netam"}}) &&
    add_test({test_model_output_shapes, "Netam Model: output shapes", {"netam"}}) &&
    add_test({test_model_batch_size_1, "Netam Model: batch size 1", {"netam"}}) &&
    add_test({test_model_rates_positive, "Netam Model: rates positive", {"netam"}}) &&
    add_test(
        {test_model_mask_zeros_output, "Netam Model: mask zeros output", {"netam"}}) &&
    add_test({test_model_wt_modifier_applied,
              "Netam Model: wt_modifier applied",
              {"netam"}}) &&
    add_test(
        {test_model_adjust_rate_bias, "Netam Model: adjust rate bias", {"netam"}}) &&
    add_test({test_model_deterministic_in_eval_mode,
              "Netam Model: deterministic in eval mode",
              {"netam"}}) &&
    add_test({test_model_different_kernel_sizes,
              "Netam Model: different kernel sizes",
              {"netam"}}) &&
    add_test({test_model_short_sequence, "Netam Model: short sequence", {"netam"}}) &&
    add_test({test_model_csp_logits_sum_behavior,
              "Netam Model: csp_logits sum behavior",
              {"netam"}}) &&
    add_test({test_model_larger_batch, "Netam Model: larger batch", {"netam"}});

#endif  // USE_NETAM

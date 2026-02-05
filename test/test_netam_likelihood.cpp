// Tests for poisson_context_log_likelihood
#include "test_common.hpp"

#ifdef USE_NETAM
#include <netam/likelihood.hpp>

#include <cmath>
#include <torch/torch.h>

using namespace netam;

namespace {

// Manually compute the log-likelihood for verification
double manual_log_likelihood(const std::vector<float>& rates,
                             const std::vector<std::vector<float>>& csp,
                             const std::vector<int64_t>& parent,
                             const std::vector<int64_t>& child) {
  // Find mutations
  std::vector<std::size_t> mut_indices;
  for (std::size_t i = 0; i < parent.size(); ++i) {
    if (parent[i] != child[i]) {
      mut_indices.push_back(i);
    }
  }

  int64_t n = static_cast<int64_t>(mut_indices.size());
  if (n == 0) {
    return 0.0;
  }

  // Sum of all rates
  double sum_rates = 0.0;
  for (float r : rates) {
    sum_rates += static_cast<double>(r);
  }

  // t_hat = n / sum_rates
  double t_hat = static_cast<double>(n) / sum_rates;

  // Sum of log(rate * csp) for mutated positions
  double log_lambda_sum = 0.0;
  for (std::size_t idx : mut_indices) {
    double rate = static_cast<double>(rates[idx]);
    double sub_prob =
        static_cast<double>(csp[idx][static_cast<std::size_t>(child[idx])]);
    log_lambda_sum += std::log(rate * sub_prob);
  }

  return log_lambda_sum + std::log(t_hat) * static_cast<double>(n) -
         static_cast<double>(n);
}

void test_no_mutations_returns_zero() {
  torch::NoGradGuard no_grad;
  // Parent and child are identical
  auto rates = torch::ones({1, 5});
  auto csp = torch::full({1, 5, 4}, 0.25f);
  auto parent = torch::tensor({0, 1, 2, 3, 0}, torch::kInt64);
  auto child = torch::tensor({0, 1, 2, 3, 0}, torch::kInt64);  // Same as parent

  auto result = poisson_context_log_likelihood(rates, csp, parent, child);

  TestAssert(result.item<double>() == 0.0);
}

void test_single_mutation() {
  torch::NoGradGuard no_grad;
  // One mutation at position 2: G(2) -> T(3)
  auto rates = torch::tensor({{1.0f, 1.0f, 2.0f, 1.0f, 1.0f}});  // [1, 5]

  // CSP: uniform 0.25 for all positions
  auto csp = torch::full({1, 5, 4}, 0.25f);

  auto parent = torch::tensor({0, 1, 2, 3, 0}, torch::kInt64);  // A C G T A
  auto child = torch::tensor({0, 1, 3, 3, 0},
                             torch::kInt64);  // A C T T A (G->T at pos 2)

  auto result = poisson_context_log_likelihood(rates, csp, parent, child);

  // Manual calculation:
  // n = 1, sum_rates = 6, t_hat = 1/6
  // rate at pos 2 = 2.0, csp[2][3] = 0.25
  // log_likelihood = log(2.0 * 0.25) + log(1/6) * 1 - 1
  //                = log(0.5) + log(1/6) - 1
  double expected = std::log(0.5) + std::log(1.0 / 6.0) - 1.0;

  TestAssert(std::abs(result.item<double>() - expected) < 1e-5);
}

void test_multiple_mutations() {
  torch::NoGradGuard no_grad;
  // Two mutations
  std::vector<float> rates_vec = {1.0f, 2.0f, 3.0f, 4.0f};
  std::vector<std::vector<float>> csp_vec = {
      {0.1f, 0.2f, 0.3f, 0.4f},
      {0.4f, 0.3f, 0.2f, 0.1f},
      {0.25f, 0.25f, 0.25f, 0.25f},
      {0.1f, 0.1f, 0.1f, 0.7f},
  };
  std::vector<int64_t> parent_vec = {0, 1, 2, 3};  // A C G T
  std::vector<int64_t> child_vec = {1, 1, 0, 3};   // C C A T (mutations at 0 and 2)

  auto rates = torch::tensor(rates_vec).unsqueeze(0);  // [1, 4]
  auto csp = torch::zeros({1, 4, 4});
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      csp[0][i][j] = csp_vec[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    }
  }
  auto parent = torch::tensor(parent_vec, torch::kInt64);
  auto child = torch::tensor(child_vec, torch::kInt64);

  auto result = poisson_context_log_likelihood(rates, csp, parent, child);
  double expected = manual_log_likelihood(rates_vec, csp_vec, parent_vec, child_vec);

  TestAssert(std::abs(result.item<double>() - expected) < 1e-5);
}

void test_all_positions_mutated() {
  torch::NoGradGuard no_grad;
  std::vector<float> rates_vec = {1.0f, 1.0f, 1.0f};
  std::vector<std::vector<float>> csp_vec = {
      {0.0f, 0.5f, 0.3f, 0.2f},  // Parent is A, can't mutate to A
      {0.3f, 0.0f, 0.4f, 0.3f},  // Parent is C
      {0.2f, 0.3f, 0.0f, 0.5f},  // Parent is G
  };
  std::vector<int64_t> parent_vec = {0, 1, 2};  // A C G
  std::vector<int64_t> child_vec = {1, 2, 3};   // C G T (all mutated)

  auto rates = torch::tensor(rates_vec).unsqueeze(0);
  auto csp = torch::zeros({1, 3, 4});
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 4; ++j) {
      csp[0][i][j] = csp_vec[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    }
  }
  auto parent = torch::tensor(parent_vec, torch::kInt64);
  auto child = torch::tensor(child_vec, torch::kInt64);

  auto result = poisson_context_log_likelihood(rates, csp, parent, child);
  double expected = manual_log_likelihood(rates_vec, csp_vec, parent_vec, child_vec);

  TestAssert(std::abs(result.item<double>() - expected) < 1e-5);
}

void test_single_position_sequence() {
  torch::NoGradGuard no_grad;
  // Edge case: sequence of length 1 with a mutation
  auto rates = torch::tensor({{2.0f}});                    // [1, 1]
  auto csp = torch::tensor({{{0.0f, 0.5f, 0.3f, 0.2f}}});  // [1, 1, 4]

  auto parent = torch::tensor({0}, torch::kInt64);  // A
  auto child = torch::tensor({1}, torch::kInt64);   // C

  auto result = poisson_context_log_likelihood(rates, csp, parent, child);

  // n = 1, sum_rates = 2, t_hat = 0.5
  // rate = 2.0, csp = 0.5
  // log_likelihood = log(2.0 * 0.5) + log(0.5) - 1 = log(1) + log(0.5) - 1
  double expected = std::log(1.0) + std::log(0.5) - 1.0;

  TestAssert(std::abs(result.item<double>() - expected) < 1e-5);
}

void test_output_is_scalar() {
  torch::NoGradGuard no_grad;
  auto rates = torch::ones({1, 10});
  auto csp = torch::full({1, 10, 4}, 0.25f);
  auto parent = torch::zeros({10}, torch::kInt64);
  auto child = torch::ones({10}, torch::kInt64);  // All mutated

  auto result = poisson_context_log_likelihood(rates, csp, parent, child);

  TestAssert(result.dim() == 0);  // Scalar tensor
}

void test_high_csp_gives_higher_likelihood() {
  torch::NoGradGuard no_grad;
  // Higher CSP for the mutation should give higher (less negative) likelihood
  auto rates = torch::ones({1, 3});
  auto parent = torch::tensor({0, 1, 2}, torch::kInt64);
  auto child = torch::tensor({1, 1, 2}, torch::kInt64);  // Mutation at pos 0: A->C

  // Low CSP for A->C
  auto csp_low = torch::full({1, 3, 4}, 0.25f);
  csp_low[0][0][1] = 0.1f;  // Low probability for A->C

  // High CSP for A->C
  auto csp_high = torch::full({1, 3, 4}, 0.25f);
  csp_high[0][0][1] = 0.9f;  // High probability for A->C

  auto result_low = poisson_context_log_likelihood(rates, csp_low, parent, child);
  auto result_high = poisson_context_log_likelihood(rates, csp_high, parent, child);

  // Higher CSP should give higher log-likelihood
  TestAssert(result_high.item<double>() > result_low.item<double>());
}

void test_high_rate_at_mutation_gives_higher_likelihood() {
  torch::NoGradGuard no_grad;
  // Higher rate at mutation site should give higher likelihood
  auto csp = torch::full({1, 3, 4}, 0.25f);
  auto parent = torch::tensor({0, 1, 2}, torch::kInt64);
  auto child = torch::tensor({1, 1, 2}, torch::kInt64);  // Mutation at pos 0

  // Low rate at mutation site
  auto rates_low = torch::ones({1, 3});
  rates_low[0][0] = 0.1f;

  // High rate at mutation site
  auto rates_high = torch::ones({1, 3});
  rates_high[0][0] = 10.0f;

  auto result_low = poisson_context_log_likelihood(rates_low, csp, parent, child);
  auto result_high = poisson_context_log_likelihood(rates_high, csp, parent, child);

  // Higher rate at mutation should give higher log-likelihood
  TestAssert(result_high.item<double>() > result_low.item<double>());
}

void test_symmetry_of_mutation_count() {
  torch::NoGradGuard no_grad;
  // Two different mutation patterns with same count should have predictable
  // relationship
  auto rates = torch::ones({1, 4});
  auto csp = torch::full({1, 4, 4}, 0.25f);

  // One mutation
  auto parent1 = torch::tensor({0, 1, 2, 3}, torch::kInt64);
  auto child1 = torch::tensor({1, 1, 2, 3}, torch::kInt64);

  // Same mutation at different position (same rates and CSPs)
  auto child2 = torch::tensor({0, 2, 2, 3}, torch::kInt64);

  auto result1 = poisson_context_log_likelihood(rates, csp, parent1, child1);
  auto result2 = poisson_context_log_likelihood(rates, csp, parent1, child2);

  // With uniform rates and CSPs, results should be equal
  TestAssert(std::abs(result1.item<double>() - result2.item<double>()) < 1e-5);
}

void test_longer_sequence() {
  torch::NoGradGuard no_grad;
  // Test with a longer sequence
  const int64_t len = 100;
  auto rates = torch::rand({1, len}) + 0.1f;  // Avoid zero rates
  auto csp = torch::rand({1, len, 4});
  // Normalize CSP to sum to 1
  csp = csp / csp.sum(-1, true);

  auto parent = torch::randint(0, 4, {len}, torch::kInt64);
  auto child = parent.clone();
  // Introduce some mutations
  child[10] = (parent[10].item<int64_t>() + 1) % 4;
  child[50] = (parent[50].item<int64_t>() + 2) % 4;
  child[90] = (parent[90].item<int64_t>() + 3) % 4;

  auto result = poisson_context_log_likelihood(rates, csp, parent, child);

  // Result should be finite
  TestAssert(std::isfinite(result.item<double>()));
}

void test_formula_components() {
  torch::NoGradGuard no_grad;
  // Verify the formula: log_lik = Σlog(λ_j) + n*log(t_hat) - n
  // where λ_j = rate_j * csp_j, t_hat = n / Σrates

  auto rates = torch::tensor({{1.0f, 2.0f, 3.0f, 4.0f}});  // sum = 10
  auto csp = torch::full({1, 4, 4}, 0.25f);

  auto parent = torch::tensor({0, 1, 2, 3}, torch::kInt64);
  auto child = torch::tensor({1, 0, 2, 3}, torch::kInt64);  // 2 mutations at pos 0,1

  auto result = poisson_context_log_likelihood(rates, csp, parent, child);

  // Manual calculation:
  // n = 2, sum_rates = 10, t_hat = 2/10 = 0.2
  // λ_0 = 1.0 * 0.25 = 0.25 (mutation A->C at pos 0)
  // λ_1 = 2.0 * 0.25 = 0.5  (mutation C->A at pos 1)
  // log_lik = log(0.25) + log(0.5) + 2*log(0.2) - 2
  double expected = std::log(0.25) + std::log(0.5) + 2.0 * std::log(0.2) - 2.0;

  TestAssert(std::abs(result.item<double>() - expected) < 1e-5);
}

}  // namespace

[[maybe_unused]] static bool reg_test_likelihood =
    add_test({test_no_mutations_returns_zero,
              "Netam Likelihood: No mutations returns zero",
              {"netam"}}) &&
    add_test({test_single_mutation, "Netam Likelihood: Single mutation", {"netam"}}) &&
    add_test(
        {test_multiple_mutations, "Netam Likelihood: Multiple mutations", {"netam"}}) &&
    add_test({test_all_positions_mutated,
              "Netam Likelihood: All positions mutated",
              {"netam"}}) &&
    add_test({test_single_position_sequence,
              "Netam Likelihood: Single position sequence",
              {"netam"}}) &&
    add_test(
        {test_output_is_scalar, "Netam Likelihood: Output is scalar", {"netam"}}) &&
    add_test({test_high_csp_gives_higher_likelihood,
              "Netam Likelihood: High CSP gives higher likelihood",
              {"netam"}}) &&
    add_test({test_high_rate_at_mutation_gives_higher_likelihood,
              "Netam Likelihood: High rate at mutation gives higher likelihood",
              {"netam"}}) &&
    add_test({test_symmetry_of_mutation_count,
              "Netam Likelihood: Symmetry of mutation count",
              {"netam"}}) &&
    add_test({test_longer_sequence, "Netam Likelihood: Longer sequence", {"netam"}}) &&
    add_test(
        {test_formula_components, "Netam Likelihood: Formula components", {"netam"}});

#endif  // USE_NETAM

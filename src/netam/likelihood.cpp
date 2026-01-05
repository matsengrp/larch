#include <netam/likelihood.hpp>

namespace netam {

torch::Tensor poisson_context_log_likelihood(
    const torch::Tensor& rates,   // [1, L]
    const torch::Tensor& csp,     // [1, L, 4]
    const torch::Tensor& parent,  // [L]
    const torch::Tensor& child    // [L]
) {
  auto rates_flat = rates.squeeze(0);  // [L]
  auto csp_flat = csp.squeeze(0);      // [L, 4]

  auto mutation_mask = parent != child;  // [L]
  int64_t n_mutations = mutation_mask.sum().item<int64_t>();

  if (n_mutations == 0) {
    return torch::tensor(0.0);
  }

  auto mutated_indices = mutation_mask.nonzero().squeeze(-1);  // [n_mutations]
  auto child_bases = child.index_select(0, mutated_indices);   // [n_mutations]

  auto rates_mut = rates_flat.index_select(0, mutated_indices);
  auto csp_mut = csp_flat.index_select(0, mutated_indices);

  auto sub_probs = csp_mut.gather(1, child_bases.unsqueeze(1)).squeeze(1);

  auto lambda_j = rates_mut * sub_probs;

  double t_hat =
      static_cast<double>(n_mutations) / rates_flat.sum().item<double>();

  return lambda_j.log().sum() +
         std::log(t_hat) * static_cast<double>(n_mutations) - n_mutations;
}

}  // namespace netam

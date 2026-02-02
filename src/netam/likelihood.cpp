#include <netam/likelihood.hpp>

namespace netam {

torch::Tensor poisson_context_log_likelihood(
    const torch::Tensor& rates,  // [1, L] or [L]
    const torch::Tensor& csp,    // [1, L, 4] or [L, 4]
    const torch::Tensor& parent, // [L']
    const torch::Tensor& child,  // [L']
    const torch::Tensor& mask    // [1, L] or [L] or empty
) {
  auto rates_flat = rates.dim() == 2 ? rates.squeeze(0) : rates;  // [L]
  auto csp_flat = csp.dim() == 3 ? csp.squeeze(0) : csp;          // [L, 4]

  // Only compare positions within the parent/child sequence length
  auto seq_len = parent.size(0);
  auto rates_seq = rates_flat.slice(0, 0, seq_len);  // [L']
  auto csp_seq = csp_flat.slice(0, 0, seq_len);      // [L', 4]

  auto mutation_mask = parent != child;  // [L']
  int64_t n_mutations = mutation_mask.sum().item<int64_t>();

  if (n_mutations == 0) {
    return torch::tensor(0.0);
  }

  auto mutated_indices = mutation_mask.nonzero().squeeze(-1);  // [n_mutations]
  auto child_bases = child.index_select(0, mutated_indices);   // [n_mutations]

  auto rates_mut = rates_seq.index_select(0, mutated_indices);
  auto csp_mut = csp_seq.index_select(0, mutated_indices);

  auto sub_probs = csp_mut.gather(1, child_bases.unsqueeze(1)).squeeze(1);

  auto lambda_j = rates_mut * sub_probs;

  // Sum rates only over valid positions (either from mask or from sequence length)
  double sum_rates;
  if (mask.defined() && mask.numel() > 0) {
    auto mask_flat = mask.dim() == 2 ? mask.squeeze(0) : mask;
    sum_rates = (rates_flat * mask_flat.to(torch::kFloat32)).sum().item<double>();
  } else {
    sum_rates = rates_seq.sum().item<double>();
  }

  double t_hat = static_cast<double>(n_mutations) / sum_rates;

  return lambda_j.log().sum() +
         std::log(t_hat) * static_cast<double>(n_mutations) - n_mutations;
}

}  // namespace netam

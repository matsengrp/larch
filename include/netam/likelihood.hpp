#pragma once

#include <torch/torch.h>

namespace netam {

// Compute log-likelihood for a parent-child pair using the Poisson context model.
// The mask indicates valid positions (1 for valid, 0 for invalid/padding).
// If mask is not provided, all positions up to the length of parent/child are used.
torch::Tensor poisson_context_log_likelihood(const torch::Tensor& rates,
                                             const torch::Tensor& csp,
                                             const torch::Tensor& parent,
                                             const torch::Tensor& child,
                                             const torch::Tensor& mask = {});

}  // namespace netam

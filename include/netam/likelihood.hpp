#pragma once

#include <torch/torch.h>

namespace netam {

torch::Tensor poisson_context_log_likelihood(const torch::Tensor& rates,
                                             const torch::Tensor& csp,
                                             const torch::Tensor& parent,
                                             const torch::Tensor& child);

}  // namespace netam

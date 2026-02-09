#pragma once

#include <torch/torch.h>
#include <netam/include-yaml.hpp>

namespace netam {

class rs_fivemer_params {
 public:
  rs_fivemer_params(std::size_t kmer_count, const YAML::Node& yaml);
  rs_fivemer_params(const rs_fivemer_params&) = default;

  std::size_t kmer_count() const noexcept;

  std::size_t kmer_length() const noexcept;

 private:
  const std::size_t kmer_count_;
  const std::size_t kmer_length_;
};

class rs_fivemer_model : public torch::nn::Module {
 public:
  rs_fivemer_model(const rs_fivemer_params& params);

  std::pair<torch::Tensor, torch::Tensor> forward(torch::Tensor encoded_parents,
                                                  torch::Tensor masks,
                                                  torch::Tensor wt_base_modifier);

  void adjust_rate_bias_by(double log_adjustment_factor);

 private:
  // R component: mutation rates per kmer
  torch::nn::Embedding r_kmer_embedding_;

  // S component: conditional substitution probabilities per kmer
  torch::nn::Embedding s_kmer_embedding_;
};

}  // namespace netam

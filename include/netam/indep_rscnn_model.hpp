#pragma once

#include <torch/torch.h>
#include <netam/include-yaml.hpp>

namespace netam {

class indep_rscnn_params {
 public:
  indep_rscnn_params(std::size_t kmer_count, const YAML::Node& yaml);
  indep_rscnn_params(const indep_rscnn_params&) = default;

  std::size_t kmer_count() const noexcept;

  std::size_t kmer_length() const noexcept;

  std::size_t embedding_dim() const noexcept;

  std::size_t filter_count() const noexcept;

  std::size_t kernel_size() const noexcept;

  double dropout_prob() const noexcept;

 private:
  const std::size_t kmer_count_;
  const std::size_t kmer_length_;
  const std::size_t embedding_dim_;
  const std::size_t filter_count_;
  const std::size_t kernel_size_;
  const double dropout_prob_;
};

class indep_rscnn_model : public torch::nn::Module {
 public:
  indep_rscnn_model(const indep_rscnn_params& params);

  std::pair<torch::Tensor, torch::Tensor> forward(torch::Tensor encoded_parents,
                                                  torch::Tensor masks,
                                                  torch::Tensor wt_base_modifier);

  void adjust_rate_bias_by(double log_adjustment_factor);

 private:
  // Asymmetric padding for "same" convolution (matches PyTorch's
  // padding="same")
  const std::size_t pad_left_;
  const std::size_t pad_right_;

  // R component layers
  torch::nn::Embedding r_kmer_embedding_;
  torch::nn::Conv1d r_conv_;
  torch::nn::Dropout r_dropout_;
  torch::nn::Linear r_linear_;

  // S component layers
  torch::nn::Embedding s_kmer_embedding_;
  torch::nn::Conv1d s_conv_;
  torch::nn::Dropout s_dropout_;
  torch::nn::Linear s_linear_;
};

}  // namespace netam

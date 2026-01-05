#include <netam/indep_rscnn_model.hpp>
#include <netam/common.hpp>

namespace netam {

/*
 * class indep_rscnn_params
 */

indep_rscnn_params::indep_rscnn_params(std::size_t kmer_count,
                                       const YAML::Node& yaml)
    : kmer_count_{kmer_count},
      kmer_length_{yaml["kmer_length"].as<std::size_t>()},
      embedding_dim_{yaml["embedding_dim"].as<std::size_t>()},
      filter_count_{yaml["filter_count"].as<std::size_t>()},
      kernel_size_{yaml["kernel_size"].as<std::size_t>()},
      dropout_prob_{yaml["dropout_prob"].as<double>()} {}

std::size_t indep_rscnn_params::kmer_count() const noexcept {
  return kmer_count_;
}

std::size_t indep_rscnn_params::kmer_length() const noexcept {
  return kmer_length_;
}

std::size_t indep_rscnn_params::embedding_dim() const noexcept {
  return embedding_dim_;
}

std::size_t indep_rscnn_params::filter_count() const noexcept {
  return filter_count_;
}

std::size_t indep_rscnn_params::kernel_size() const noexcept {
  return kernel_size_;
}

double indep_rscnn_params::dropout_prob() const noexcept {
  return dropout_prob_;
}

/*
 * class indep_rscnn_model
 */

indep_rscnn_model::indep_rscnn_model(const indep_rscnn_params& params)
    :  // Calculate asymmetric padding for "same" convolution (matches PyTorch)
       // total_padding = kernel_size - 1
       // pad_left = total_padding // 2
       // pad_right = total_padding - pad_left
      pad_left_{(params.kernel_size() - 1) / 2},
      pad_right_{params.kernel_size() - 1 - pad_left_},

      // R component layers (no padding in conv, we pad manually)
      r_kmer_embedding_{register_module(
          "r_kmer_embedding",
          torch::nn::Embedding{params.kmer_count(), params.embedding_dim()})},
      r_conv_{
          register_module("r_conv", torch::nn::Conv1d{torch::nn::Conv1dOptions{
                                        signed_cast(params.embedding_dim()),
                                        signed_cast(params.filter_count()),
                                        signed_cast(params.kernel_size())}})},
      r_dropout_{register_module("r_dropout",
                                 torch::nn::Dropout{params.dropout_prob()})},
      r_linear_{register_module("r_linear",
                                torch::nn::Linear{params.filter_count(), 1})},

      // S component layers (no padding in conv, we pad manually)
      s_kmer_embedding_{register_module(
          "s_kmer_embedding",
          torch::nn::Embedding{params.kmer_count(), params.embedding_dim()})},
      s_conv_{
          register_module("s_conv", torch::nn::Conv1d{torch::nn::Conv1dOptions(
                                        signed_cast(params.embedding_dim()),
                                        signed_cast(params.filter_count()),
                                        signed_cast(params.kernel_size()))})},
      s_dropout_{register_module("s_dropout",
                                 torch::nn::Dropout{params.dropout_prob()})},
      s_linear_{register_module("s_linear",
                                torch::nn::Linear{params.filter_count(), 4})} {}

std::pair<torch::Tensor, torch::Tensor> indep_rscnn_model::forward(
    torch::Tensor encoded_parents, torch::Tensor masks,
    torch::Tensor wt_base_modifier) {
  // Process R component
  auto r_kmer_embeds = r_kmer_embedding_->forward(encoded_parents);
  r_kmer_embeds = r_kmer_embeds.permute({0, 2, 1});  // [B, E, L]
  // Apply asymmetric "same" padding before convolution
  r_kmer_embeds = torch::nn::functional::pad(
      r_kmer_embeds, torch::nn::functional::PadFuncOptions(
                         {signed_cast(pad_left_), signed_cast(pad_right_)}));
  auto r_conv_out = torch::relu(r_conv_->forward(r_kmer_embeds));
  r_conv_out = r_dropout_->forward(r_conv_out);
  r_conv_out = r_conv_out.permute({0, 2, 1});  // [B, L, F]

  auto log_rates = r_linear_->forward(r_conv_out).squeeze(-1);  // [B, L]
  auto rates = torch::exp(log_rates * masks);

  // Process S component
  auto s_kmer_embeds = s_kmer_embedding_->forward(encoded_parents);
  s_kmer_embeds = s_kmer_embeds.permute({0, 2, 1});  // [B, E, L]
  // Apply asymmetric "same" padding before convolution
  s_kmer_embeds = torch::nn::functional::pad(
      s_kmer_embeds, torch::nn::functional::PadFuncOptions(
                         {signed_cast(pad_left_), signed_cast(pad_right_)}));
  auto s_conv_out = torch::relu(s_conv_->forward(s_kmer_embeds));
  s_conv_out = s_dropout_->forward(s_conv_out);
  s_conv_out = s_conv_out.permute({0, 2, 1});  // [B, L, F]

  auto csp_logits = s_linear_->forward(s_conv_out);  // [B, L, 4]
  csp_logits = csp_logits * masks.unsqueeze(-1);
  csp_logits = csp_logits + wt_base_modifier;

  return {rates, csp_logits};
}

void indep_rscnn_model::adjust_rate_bias_by(double log_adjustment_factor) {
  r_linear_->bias.data() += log_adjustment_factor;
}

}  // namespace netam

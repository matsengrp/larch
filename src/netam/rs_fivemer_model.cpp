#include <netam/rs_fivemer_model.hpp>
#include <netam/common.hpp>

namespace netam {

/*
 * class rs_fivemer_params
 */

rs_fivemer_params::rs_fivemer_params(std::size_t kmer_count,
                                     const YAML::Node& yaml)
    : kmer_count_{kmer_count},
      kmer_length_{yaml["kmer_length"].as<std::size_t>()} {
  Assert(kmer_length_ == 5);
}

std::size_t rs_fivemer_params::kmer_count() const noexcept {
  return kmer_count_;
}

std::size_t rs_fivemer_params::kmer_length() const noexcept {
  return kmer_length_;
}

/*
 * class rs_fivemer_model
 */

rs_fivemer_model::rs_fivemer_model(const rs_fivemer_params& params)
    : r_kmer_embedding_{register_module(
          "r_kmer_embedding",
          torch::nn::Embedding{params.kmer_count(), 1})},
      s_kmer_embedding_{register_module(
          "s_kmer_embedding",
          torch::nn::Embedding{params.kmer_count(), 4})} {}

std::pair<torch::Tensor, torch::Tensor> rs_fivemer_model::forward(
    torch::Tensor encoded_parents, torch::Tensor masks,
    torch::Tensor wt_base_modifier) {
  // R component: get log rates from embedding and exponentiate
  auto log_kmer_rates = r_kmer_embedding_->forward(encoded_parents).squeeze(-1);
  auto rates = torch::exp(log_kmer_rates * masks);

  // S component: get CSP logits from embedding
  auto csp_logits = s_kmer_embedding_->forward(encoded_parents);

  // When we have an N, set all the CSP logits to 0, resulting in a uniform
  // prediction
  csp_logits = csp_logits * masks.unsqueeze(-1);

  // Apply wt_base_modifier to zero out wild-type base predictions
  csp_logits = csp_logits + wt_base_modifier;

  return {rates, csp_logits};
}

void rs_fivemer_model::adjust_rate_bias_by(double log_adjustment_factor) {
  torch::NoGradGuard no_grad;
  r_kmer_embedding_->weight.data() += log_adjustment_factor;
}

}  // namespace netam

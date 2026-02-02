#pragma once

#include <netam/kmer_sequence_encoder.hpp>
#include <netam/indep_rscnn_model.hpp>
#include <netam/rs_fivemer_model.hpp>

#include <filesystem>
#include <variant>

namespace netam {

// Unified model interface that wraps the variant
class model_wrapper {
 public:
  using model_variant = std::variant<indep_rscnn_model, rs_fivemer_model>;

  model_wrapper(model_variant model);

  std::pair<torch::Tensor, torch::Tensor> forward(
      torch::Tensor encoded_parents, torch::Tensor masks,
      torch::Tensor wt_base_modifier);

  void eval();

  void adjust_rate_bias_by(double log_adjustment_factor);

  // Allow access for state dict loading
  model_variant& get_model() noexcept { return model_; }

 private:
  model_variant model_;
};

class crepe {
 public:
  crepe(const std::filesystem::path& weights_path,
        const std::filesystem::path& yaml_path);

  kmer_sequence_encoder& encoder() noexcept;

  model_wrapper* operator->() noexcept;

 private:
  YAML::Node yaml_;
  kmer_sequence_encoder encoder_;
  model_wrapper model_;
};

}  // namespace netam

#include <netam/crepe.hpp>
#include <netam/common.hpp>

#include <fstream>

namespace netam {

/*
 * class model_wrapper
 */

model_wrapper::model_wrapper(model_variant model) : model_{std::move(model)} {}

std::pair<torch::Tensor, torch::Tensor> model_wrapper::forward(
    torch::Tensor encoded_parents, torch::Tensor masks,
    torch::Tensor wt_base_modifier) {
  return std::visit(
      [&](auto& m) { return m.forward(encoded_parents, masks, wt_base_modifier); },
      model_);
}

void model_wrapper::eval() {
  std::visit([](auto& m) { m.eval(); }, model_);
}

void model_wrapper::adjust_rate_bias_by(double log_adjustment_factor) {
  std::visit([&](auto& m) { m.adjust_rate_bias_by(log_adjustment_factor); }, model_);
}

/*
 * class crepe
 */

namespace {

model_wrapper::model_variant create_model(std::size_t kmer_count,
                                          const YAML::Node& yaml) {
  std::string model_class = yaml["model_class"].as<std::string>();

  if (model_class == "IndepRSCNNModel") {
    return indep_rscnn_model{
        indep_rscnn_params{kmer_count, yaml["model_hyperparameters"]}};
  } else if (model_class == "RSFivemerModel") {
    return rs_fivemer_model{
        rs_fivemer_params{kmer_count, yaml["model_hyperparameters"]}};
  } else {
    throw std::runtime_error("Unknown model_class: " + model_class);
  }
}

void load_state_dict(model_wrapper::model_variant& model,
                     const std::filesystem::path& weights_path) {
  std::ifstream file(weights_path, std::ios::binary);
  std::vector<char> data{std::istreambuf_iterator<char>{file},
                         std::istreambuf_iterator<char>{}};
  auto loaded = torch::pickle_load(data);
  auto state_dict = loaded.toGenericDict();
  torch::NoGradGuard no_grad;

  std::visit(
      [&state_dict](auto& m) {
        for (auto& param : m.named_parameters()) {
          std::string name = param.key();
          Assert(state_dict.contains(name));
          torch::Tensor loaded_tensor = state_dict.at(name).toTensor();
          param.value().copy_(loaded_tensor);
        }

        for (auto& buffer : m.named_buffers()) {
          std::string name = buffer.key();
          Assert(state_dict.contains(name));
          torch::Tensor loaded_tensor = state_dict.at(name).toTensor();
          buffer.value().copy_(loaded_tensor);
        }
      },
      model);
}

}  // namespace

crepe::crepe(const std::filesystem::path& weights_path,
             const std::filesystem::path& yaml_path)
    : yaml_{YAML::LoadFile(yaml_path)},
      encoder_{yaml_["encoder_parameters"]},
      model_{create_model(encoder_.kmer_count(), yaml_)} {
  load_state_dict(model_.get_model(), weights_path);
}

kmer_sequence_encoder& crepe::encoder() noexcept { return encoder_; }

model_wrapper* crepe::operator->() noexcept { return &model_; }

}  // namespace netam

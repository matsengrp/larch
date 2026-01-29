#include <netam/crepe.hpp>
#include <netam/common.hpp>

#include <fstream>

namespace netam {

crepe::crepe(const std::filesystem::path& weights_path,
             const std::filesystem::path& yaml_path)
    : yaml_{YAML::LoadFile(yaml_path)},
      encoder_{yaml_["encoder_parameters"]},
      model_{{encoder_.kmer_count(), yaml_["model_hyperparameters"]}} {
  std::ifstream file(weights_path, std::ios::binary);
  std::vector<char> data{std::istreambuf_iterator<char>{file},
                         std::istreambuf_iterator<char>{}};
  auto loaded = torch::pickle_load(data);
  auto state_dict = loaded.toGenericDict();
  torch::NoGradGuard no_grad;

  for (auto& param : model_.named_parameters()) {
    std::string name = param.key();
    Assert(state_dict.contains(name));
    torch::Tensor loaded_tensor = state_dict.at(name).toTensor();
    param.value().copy_(loaded_tensor);
  }

  for (auto& buffer : model_.named_buffers()) {
    std::string name = buffer.key();
    Assert(state_dict.contains(name));
    torch::Tensor loaded_tensor = state_dict.at(name).toTensor();
    buffer.value().copy_(loaded_tensor);
  }
}

kmer_sequence_encoder& crepe::encoder() noexcept { return encoder_; }

indep_rscnn_model* crepe::operator->() noexcept { return &model_; }

}  // namespace netam

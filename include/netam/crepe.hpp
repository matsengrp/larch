#pragma once

#include <netam/kmer_sequence_encoder.hpp>
#include <netam/indep_rscnn_model.hpp>

#include <filesystem>

namespace netam {

class crepe {
 public:
  crepe(const std::filesystem::path& weights_path,
        const std::filesystem::path& yaml_path);

  kmer_sequence_encoder& encoder() noexcept;

  indep_rscnn_model* operator->() noexcept;

 private:
  YAML::Node yaml_;
  kmer_sequence_encoder encoder_;
  indep_rscnn_model model_;
};

}  // namespace netam

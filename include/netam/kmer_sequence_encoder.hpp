#pragma once

#include <torch/torch.h>
#include <netam/include-yaml.hpp>

#include <vector>
#include <string>
#include <unordered_map>

namespace netam {

class kmer_sequence_encoder {
 public:
  kmer_sequence_encoder(const YAML::Node& yaml);

  std::pair<torch::Tensor, torch::Tensor> encode_sequence(
      const std::string& sequence) const;

  std::size_t kmer_count() const noexcept;

  std::size_t kmer_length() const noexcept;

  std::size_t site_count() const noexcept;

  // Encode a sequence as base indices (A=0, C=1, G=2, T=3, other=4)
  // This is used for likelihood calculations, not for model input.
  static torch::Tensor encode_bases(const std::string& sequence);

 private:
  const std::size_t kmer_length_;
  const std::size_t site_count_;
  const std::size_t overhang_length_;
  const std::vector<std::string> all_kmers_;
  const std::unordered_map<std::string, std::size_t> kmer_to_index_;

  static constexpr const char* BASES = "ACGT";
  static constexpr float BIG = 1e9f;

  static std::vector<std::string> generate_kmers(std::size_t length);

  torch::Tensor compute_wt_base_modifier(const std::string& parent) const;

  static std::size_t get_base_index(char base);

  static std::string to_upper(const std::string& str);
};

}  // namespace netam

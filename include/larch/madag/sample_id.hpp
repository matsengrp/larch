#pragma once

struct SampleId {
  std::optional<std::string> sample_id_;

  inline bool empty() const {
    if (not sample_id_.has_value()) {
      return false;
    }
    return sample_id_.value().empty();
  }

  inline std::string ToString() const { return sample_id_.value_or(std::string{}); }
  inline size_t Hash() const {
    return std::hash<std::string>{}(sample_id_.value_or(std::string{}));
  }
  inline static const SampleId* Empty();
};

template <>
struct std::hash<SampleId> {
  inline std::size_t operator()(const SampleId& sid) const noexcept;
};

template <>
struct std::equal_to<SampleId> {
  inline bool operator()(const SampleId& lhs, const SampleId& rhs) const noexcept;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<SampleId, CRTP, Tag> {
  const std::optional<std::string>& GetSampleId() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<SampleId, CRTP, Tag> {
  void SetSampleId(const std::optional<std::string>& sample_id) const;
};

#include "larch/impl/madag/sample_id_impl.hpp"

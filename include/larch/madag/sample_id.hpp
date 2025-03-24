#pragma once

struct SampleId {
  MOVE_ONLY(SampleId);

  SampleId() = default;
  SampleId(std::optional<std::string> id)
      : hash_{[](const std::optional<std::string>& x) -> size_t {
          if (not x.has_value()) {
            return 0;
          } else {
            return std::hash<std::string>{}(x.value());
          }
        }(id)},
        sample_id_{std::move(id)} {}

  size_t hash_ = 0;
  std::optional<std::string> sample_id_;

  template <typename CRTP>
  inline SampleId Copy(const CRTP*) const {
    return SampleId(sample_id_);
  }

  inline bool empty() const {
    if (not sample_id_.has_value()) {
      return true;
    }
    return sample_id_.value().empty();
  }

  inline std::string ToString() const { return sample_id_.value_or(std::string{}); }
  inline size_t Hash() const { return hash_; }
  inline static const SampleId* GetEmpty();
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

  inline bool HaveSampleId() const { return not GetFeatureStorage(this).get().empty(); }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<SampleId, CRTP, Tag> {
  void SetSampleId(const std::optional<std::string>& sample_id) const;
};

#include "larch/impl/madag/sample_id_impl.hpp"

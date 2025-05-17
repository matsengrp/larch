#pragma once

class SampleIdStorage;

template <>
struct std::hash<SampleIdStorage> {
  inline std::size_t operator()(const SampleIdStorage& sid) const noexcept;
};

template <>
struct std::equal_to<SampleIdStorage> {
  inline bool operator()(const SampleIdStorage& lhs,
                         const SampleIdStorage& rhs) const noexcept;
};

class SampleIdStorage {
 public:
  MOVE_ONLY(SampleIdStorage);

  inline size_t Hash() const { return hash_; }
  inline const std::string& Value() const { return value_; }

  static const SampleIdStorage& Get(std::string&& x) {
    return *values_.emplace(SampleIdStorage{std::move(x)}).first;
  }

 private:
  SampleIdStorage(std::string&& value)
      : hash_{std::hash<std::string>{}(value)}, value_{std::move(value)} {}

  static inline std::unordered_set<SampleIdStorage> values_;

  const size_t hash_;
  const std::string value_;
};

struct SampleId {
  MOVE_ONLY(SampleId);

  SampleId() : target_{nullptr} {}

  static SampleId Make(std::string_view id) { return std::string{id}; }

  template <typename CRTP>
  SampleId Copy(const CRTP*) const {
    return SampleId(target_);
  }

  inline bool empty() const {
    if (not target_) {
      return true;
    }
    return target_->Value().empty();
  }

  inline std::string ToString() const {
    if (not target_) {
      return {};
    }
    return std::string{target_->Value()};
  }

  inline size_t Hash() const {
    if (not target_) {
      return 0;
    }
    return target_->Hash();
  }

  inline static const SampleId* GetEmpty();

 private:
  template <typename>
  friend struct std::equal_to;
  template <typename, typename, typename>
  friend struct FeatureConstView;
  template <typename, typename, typename>
  friend struct FeatureMutableView;

  SampleId(std::string&& x) : target_{&SampleIdStorage::Get(std::move(x))} {}
  SampleId(const SampleIdStorage* x) : target_{x} {}

  const SampleIdStorage* target_;
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
  std::optional<std::string_view> GetSampleId() const;

  inline bool HaveSampleId() const { return not GetFeatureStorage(this).get().empty(); }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<SampleId, CRTP, Tag> {
  void SetSampleId(std::optional<std::string_view> sample_id) const;
};

#include "larch/impl/madag/sample_id_impl.hpp"

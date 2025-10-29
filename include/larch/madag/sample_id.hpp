#pragma once

#include <boost/unordered/unordered_set.hpp>

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

/**
 * @brief Internal storage for unique sample identifiers with efficient hash-based
 * lookup.
 *
 * SampleIdStorage implements a flyweight pattern for sample IDs, ensuring that each
 * unique sample identifier string is stored only once in memory. It provides fast
 * hash-based comparison and lookup operations. This class is used internally by
 * SampleId and should not be used directly.
 */
class SampleIdStorage {
 public:
  inline size_t Hash() const { return hash_; }
  inline const std::string& Value() const { return value_; }

  static const SampleIdStorage& Get(std::string_view x) {
    std::shared_lock read_lock{mtx_};
    auto result = values_.find(x);
    if (result != values_.end()) {
      return *result;
    }
    read_lock.unlock();
    std::unique_lock write_lock{mtx_};
    return *values_.emplace(SampleIdStorage{std::string{x}}).first;
  }

 private:
  SampleIdStorage(std::string value)
      : hash_{std::hash<std::string>{}(value)}, value_{std::move(value)} {}

  struct key_hash {
    using is_transparent = void;
    size_t operator()(const SampleIdStorage& x) const { return x.Hash(); };
    size_t operator()(std::string_view x) const {
      return std::hash<std::string_view>{}(x);
    };
  };

  struct key_equal {
    using is_transparent = void;
    bool operator()(const SampleIdStorage& lhs, const SampleIdStorage& rhs) const {
      return lhs.Value() == rhs.Value();
    }

    bool operator()(std::string_view lhs, const SampleIdStorage& rhs) const {
      return lhs == rhs.Value();
    }
  };

  static inline boost::unordered_set<SampleIdStorage, key_hash, key_equal> values_{
      1000, key_hash{}, key_equal{}};
  static inline std::shared_mutex mtx_;

  const size_t hash_;
  const std::string value_;
};

struct SampleId;
using UniqueData = SampleId;

/**
 * @brief Lightweight identifier for samples in phylogenetic trees.
 *
 * SampleId provides an efficient way to represent and compare sample identifiers in
 * phylogenetic data structures. It uses the flyweight pattern through SampleIdStorage
 * to ensure memory efficiency when dealing with many samples. The class supports fast
 * equality comparison through pointer comparison and provides hash support for use in
 * standard containers.
 */
struct SampleId {
  SampleId(SampleId&&) noexcept = default;
  SampleId(const SampleId&) noexcept = default;
  SampleId& operator=(SampleId&&) noexcept = default;
  SampleId& operator=(const SampleId&) noexcept = default;

  SampleId() : target_{nullptr} {}

  static SampleId Make(std::string_view id) { return id; }

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

  inline static UniqueData GetEmpty();

 private:
  template <typename>
  friend struct std::equal_to;
  template <typename>
  friend struct std::less;
  template <typename, typename, typename>
  friend struct FeatureConstView;
  template <typename, typename, typename>
  friend struct FeatureMutableView;

  friend bool operator==(const SampleId& lhs, const SampleId& rhs) noexcept;
  friend bool operator!=(const SampleId& lhs, const SampleId& rhs) noexcept;
  friend bool operator<(const SampleId& lhs, const SampleId& rhs) noexcept;
  friend bool operator>(const SampleId& lhs, const SampleId& rhs) noexcept;
  friend bool operator<=(const SampleId& lhs, const SampleId& rhs) noexcept;
  friend bool operator>=(const SampleId& lhs, const SampleId& rhs) noexcept;

  SampleId(std::string_view x) : target_{&SampleIdStorage::Get(x)} {}
  SampleId(const SampleIdStorage* x) : target_{x} {}

  const SampleIdStorage* target_;
};

inline bool operator==(const SampleId& lhs, const SampleId& rhs) noexcept {
  return lhs.target_ == rhs.target_;
}

inline bool operator!=(const SampleId& lhs, const SampleId& rhs) noexcept {
  return lhs.target_ != rhs.target_;
}

inline bool operator<(const SampleId& lhs, const SampleId& rhs) noexcept {
  return lhs.target_ < rhs.target_;
}

inline bool operator>(const SampleId& lhs, const SampleId& rhs) noexcept {
  return lhs.target_ > rhs.target_;
}

inline bool operator<=(const SampleId& lhs, const SampleId& rhs) noexcept {
  return lhs.target_ <= rhs.target_;
}

inline bool operator>=(const SampleId& lhs, const SampleId& rhs) noexcept {
  return lhs.target_ >= rhs.target_;
}

template <>
struct std::hash<SampleId> {
  inline std::size_t operator()(const SampleId& sid) const noexcept;
};

template <>
struct std::equal_to<SampleId> {
  inline bool operator()(const SampleId& lhs, const SampleId& rhs) const noexcept;
};

template <>
struct std::less<SampleId> {
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

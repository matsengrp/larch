#pragma once

#include <functional>
#include <memory>
#include <mutex>
#include <unordered_set>
#include <optional>

#include "larch/common.hpp"

template <typename T>
class ConcurrentUnorderedSet {
 public:
  ConcurrentUnorderedSet() = default;

  ConcurrentUnorderedSet(const ConcurrentUnorderedSet& other) {
    std::unique_lock lock{other.mtx_};
    data_ = other.data_;
  }

  ConcurrentUnorderedSet(ConcurrentUnorderedSet&& other) {
    std::unique_lock lock{other.mtx_};
    data_ = std::move(other.data_);
  }

  ConcurrentUnorderedSet& operator=(const ConcurrentUnorderedSet& other) {
    std::unique_lock other_lock{other.mtx_};
    std::unique_lock lock{mtx_};
    data_ = other.data_;
    return *this;
  }

  ConcurrentUnorderedSet& operator=(ConcurrentUnorderedSet&& other) {
    std::unique_lock other_lock{other.mtx_};
    std::unique_lock lock{mtx_};
    data_ = std::move(other.data_);
    return *this;
  }

  bool Contains(const T& value) const {
    std::unique_lock lock{mtx_};
    return data_.find(value) != data_.end();
  }

  std::optional<std::reference_wrapper<const T>> Find(const T& value) {
    std::unique_lock lock{mtx_};
    auto result = data_.find(value);
    if (result == data_.end()) {
      return std::nullopt;
    }
    return std::ref(*result);
  }

  std::pair<std::reference_wrapper<const T>, bool> Insert(T&& value) {
    std::unique_lock lock{mtx_};
    auto result = data_.insert(std::forward<T>(value));
    return {std::ref(*result.first), result.second};
  }

 private:
  std::unordered_set<T> data_;
  mutable std::mutex mtx_;
};

#include "larch/impl/parallel/node_hashset_impl.hpp"

#pragma once

#include <functional>
#include <memory>
#include <unordered_set>
#include <optional>

#include "larch/parallel/parallel_common.hpp"

template <typename T, typename M = std::shared_mutex>
class ConcurrentUnorderedSet {
 public:
  ConcurrentUnorderedSet() = default;

  ConcurrentUnorderedSet(const ConcurrentUnorderedSet& other) {
    SharedLock<M> lock{other.mtx_};
    data_ = other.data_;
  }

  ConcurrentUnorderedSet(ConcurrentUnorderedSet&& other) {
    SharedLock<M> lock{other.mtx_};
    data_ = std::move(other.data_);
  }

  ConcurrentUnorderedSet& operator=(const ConcurrentUnorderedSet& other) {
    SharedLock<M> other_lock{other.mtx_};
    std::unique_lock lock{mtx_};
    data_ = other.data_;
    return *this;
  }

  ConcurrentUnorderedSet& operator=(ConcurrentUnorderedSet&& other) {
    SharedLock<M> other_lock{other.mtx_};
    std::unique_lock lock{mtx_};
    data_ = std::move(other.data_);
    return *this;
  }

  bool Contains(const T& value) const {
    SharedLock<M> lock{mtx_};
    return data_.find(value) != data_.end();
  }

  std::optional<Accessor<const T, M>> Find(const T& value) {
    SharedLock<M> lock{mtx_};
    auto result = data_.find(value);
    if (result == data_.end()) {
      return std::nullopt;
    }
    return Accessor<const T, M>{*result, mtx_};
  }

  std::pair<Accessor<const T, M>, bool> Insert(T&& value) {
    std::unique_lock lock{mtx_};
    auto result = data_.insert(std::forward<T>(value));
    return {Accessor<const T, M>{*result.first, mtx_}, result.second};
  }

 private:
  std::unordered_set<T> data_;
  mutable M mtx_;
};

#include "larch/impl/parallel/node_hashset_impl.hpp"

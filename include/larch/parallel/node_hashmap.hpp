#pragma once

#include <functional>
#include <memory>
#include <unordered_map>
#include <optional>

#include "larch/parallel/parallel_common.hpp"

template <typename K, typename V, typename M = std::shared_mutex>
class ConcurrentUnorderedMap {
 public:
  ConcurrentUnorderedMap() = default;

  ConcurrentUnorderedMap(const ConcurrentUnorderedMap& other) {
    SharedLock<M> lock{other.mtx_};
    data_ = other.data_;
  }

  ConcurrentUnorderedMap(ConcurrentUnorderedMap&& other) {
    SharedLock<M> lock{other.mtx_};
    data_ = std::move(other.data_);
  }

  ConcurrentUnorderedMap& operator=(const ConcurrentUnorderedMap& other) {
    SharedLock<M> other_lock{other.mtx_};
    std::unique_lock lock{mtx_};
    data_ = other.data_;
    return *this;
  }

  ConcurrentUnorderedMap& operator=(ConcurrentUnorderedMap&& other) {
    SharedLock<M> other_lock{other.mtx_};
    std::unique_lock lock{mtx_};
    data_ = std::move(other.data_);
    return *this;
  }

  using Value = std::pair<K, V>;

  size_t Size() const {
    SharedLock<M> lock{mtx_};
    return data_.size();
  }

  Accessor<V, M> AtDefault(const K& key) {
    std::unique_lock lock{mtx_};
    return {data_[key], mtx_};
  }

  Accessor<V, M> At(const K& key) {
    std::unique_lock lock{mtx_};
    return {data_.at(key), mtx_};
  }

  Accessor<const V, M> At(const K& key) const {
    SharedLock<M> lock{mtx_};
    return {data_.at(key), mtx_};
  }

  auto All() {
    std::unique_lock lock{mtx_};
    return data_ | ranges::views::all;
  }

  auto All() const {
    SharedLock<M> lock{mtx_};
    return data_ | ranges::views::all;
  }

  std::optional<Accessor<V, M>> Find(const K& key) {
    std::unique_lock lock{mtx_};
    auto result = data_.find(key);
    if (result == data_.end()) {
      return std::nullopt;
    }
    return Accessor<V, M>{result->second, mtx_};
  }

  std::pair<Accessor<V, M>, bool> Insert(Value&& value) {
    std::unique_lock lock{mtx_};
    auto result = data_.insert(std::forward<Value>(value));
    return {{result.first->second, mtx_}, result.second};
  }

  template <typename... Args>
  std::pair<Accessor<V, M>, bool> Emplace(Args&&... args) {
    std::unique_lock lock{mtx_};
    auto result = data_.emplace(std::forward<Args>(args)...);
    return {{result.first->second, mtx_}, result.second};
  }

  void Clear() {
    std::unique_lock lock{mtx_};
    data_.clear();
  }

 private:
  std::unordered_map<K, V> data_;
  mutable M mtx_;
};

#include "larch/impl/parallel/node_hashmap_impl.hpp"

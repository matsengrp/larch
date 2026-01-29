#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "larch/fixed_array.hpp"
#include "larch/parallel/parallel_common.hpp"

template <typename K, typename V>
class GrowableHashMap {
  using Container = std::unordered_map<K, V>;

 public:
  using value_type = std::pair<K, V>;

  GrowableHashMap(GrowableHashMap&& other) = default;
  GrowableHashMap& operator=(GrowableHashMap&& other) = default;

  explicit GrowableHashMap(size_t buckets) : buckets_{buckets} {}

  std::pair<V&, bool> insert(value_type&& value) {
    size_t hash = std::hash<K>{}(value.first);
    auto& [mutex, data] = buckets_.at(hash % buckets_.size());

    auto wlock = WriteLock(mutex);
    auto result = data.insert(std::forward<value_type>(value));
    return {result.first->second, result.second};
  }

  template <typename Fn>
  void insert(value_type&& value, Fn&& fn) {
    size_t hash = std::hash<K>{}(value.first);
    auto& [mutex, data] = buckets_.at(hash % buckets_.size());

    auto wlock = WriteLock(mutex);
    auto result = data.insert(std::forward<value_type>(value));
    std::invoke(std::forward<Fn>(fn),
                std::pair<V&, bool>{result.first->second, result.second});
  }

  template <typename T>
  std::pair<V&, bool> insert_or_assign(K&& k, T&& v) {
    size_t hash = std::hash<K>{}(k);
    auto& [mutex, data] = buckets_.at(hash % buckets_.size());

    auto wlock = WriteLock(mutex);
    auto result = data.insert_or_assign(std::forward<K>(k), std::forward<T>(v));
    return {result.first->second, result.second};
  }

  std::pair<V&, bool> insert_or_assign(const K& k, const V& v) {
    size_t hash = std::hash<K>{}(k);
    auto& [mutex, data] = buckets_.at(hash % buckets_.size());

    auto wlock = WriteLock(mutex);
    auto result = data.insert_or_assign(k, v);
    return {result.first->second, result.second};
  }

  const V* find(const K& k) const {
    size_t hash = std::hash<K>{}(k);
    auto& [mutex, data] = buckets_.at(hash % buckets_.size());

    auto rlock = ReadLock(mutex);
    auto result = data.find(k);
    if (result == data.end()) {
      return nullptr;
    }
    return std::addressof(result->second);
  }

  V* find(const K& k) {
    size_t hash = std::hash<K>{}(k);
    auto& [mutex, data] = buckets_.at(hash % buckets_.size());

    auto rlock = ReadLock(mutex);
    auto result = data.find(k);
    if (result == data.end()) {
      return nullptr;
    }
    return std::addressof(result->second);
  }

  const V& at(const K& k) const {
    auto* result = find(k);
    Assert(result != nullptr);
    return *result;
  }

  V& at(const K& k) {
    auto* result = find(k);
    Assert(result != nullptr);
    return *result;
  }

  template <typename F, typename... Args>
  decltype(auto) ReadAll(F&& func, Args&&... args) const {
    // TODO move to member
    std::vector<std::shared_lock<std::shared_mutex>> locks;
    locks.reserve(buckets_.size());
    for (auto& i : buckets_) {
      locks.push_back(ReadLock(i.mutex_));
    }
    auto data =
        buckets_ |
        ranges::views::transform([](const auto& i) -> const auto& { return i.data_; }) |
        ranges::views::join;
    if constexpr (std::is_void_v<decltype(std::invoke(std::forward<F>(func), data,
                                                      std::forward<Args>(args)...))>) {
      std::invoke(std::forward<F>(func), data, std::forward<Args>(args)...);
    } else {
      return std::invoke(std::forward<F>(func), data, std::forward<Args>(args)...);
    }
  }

  [[nodiscard]] size_t size() const {
    // TODO move to member
    std::vector<std::shared_lock<std::shared_mutex>> locks;
    locks.reserve(buckets_.size());
    for (auto& i : buckets_) {
      locks.push_back(ReadLock(i.mutex_));
    }
    size_t result = 0;
    for (auto& i : buckets_) {
      result += i.data_.size();
    }
    return result;
  }

 private:
  struct Bucket {
    mutable std::shared_mutex mutex_;
    Container data_;
  };

  FixedArray<Bucket> buckets_;
};

template <typename K>
class GrowableHashSet {
 public:
  GrowableHashSet(GrowableHashSet&& other) = default;
  GrowableHashSet& operator=(GrowableHashSet&& other) = default;

  explicit GrowableHashSet(size_t buckets) : buckets_{buckets} {}

  template <typename T>
  std::pair<const K&, bool> insert(T&& value) {
    size_t hash = std::hash<K>{}(value);
    auto& [mutex, data] = buckets_.at(hash % buckets_.size());

    auto wlock = WriteLock(mutex);
    auto result = data.insert(std::forward<T>(value));
    return {*result.first, result.second};
  }

  template <typename Fn>
  std::pair<const K&, bool> insert(const K& value, Fn&& maker) {
    size_t hash = std::hash<K>{}(value);
    auto& [mutex, data] = buckets_.at(hash % buckets_.size());

    auto rlock = ReadLock(mutex);
    auto it = data.find(value);
    if (it == data.end()) {
      rlock.unlock();
      auto wlock = WriteLock(mutex);  // TODO upgradeable lock
      size_t size = data.size();
      auto result =
          data.insert(data.end(), std::invoke(std::forward<Fn>(maker), value));
      return {*result, data.size() != size};
    } else {
      return {*it, false};
    }
  }

  const K* find(const K& k) const {
    size_t hash = std::hash<K>{}(k);
    auto& [mutex, data] = buckets_.at(hash % buckets_.size());

    auto rlock = ReadLock(mutex);
    auto result = data.find(k);
    if (result == data.end()) {
      return nullptr;
    }
    return std::addressof(*result);
  }

 private:
  struct Bucket {
    mutable std::shared_mutex mutex_;
    std::unordered_set<K> data_;
  };

  FixedArray<Bucket> buckets_;
};

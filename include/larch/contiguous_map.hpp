#pragma once

#include <vector>
#include <algorithm>

template <typename K, typename V>
class ContiguousMap {
 public:
  ContiguousMap() = default;
  ContiguousMap(ContiguousMap&&) = default;
  ContiguousMap& operator=(ContiguousMap&&) = default;
  ContiguousMap& operator=(const ContiguousMap&) = delete;

  ContiguousMap Copy() const { return ContiguousMap{*this}; }

  auto begin() const { return data_.begin(); }

  auto end() const { return data_.end(); }

  auto LowerBound(K key) const {
    return std::lower_bound(data_.begin(), data_.end(), key,
                            [](std::pair<K, V> lhs, K rhs) { return lhs.first < rhs; });
  }

  bool operator==(const ContiguousMap& other) const { return data_ == other.data_; }

  bool operator!=(const ContiguousMap& other) const { return data_ != other.data_; }

  bool operator<(const ContiguousMap& other) const { return data_ < other.data_; }

  bool empty() const { return data_.empty(); }

  size_t size() const { return data_.size(); }

  auto begin() { return data_.begin(); }

  auto end() { return data_.end(); }

  auto LowerBound(K key) {
    return std::lower_bound(data_.begin(), data_.end(), key,
                            [](std::pair<K, V> lhs, K rhs) { return lhs.first < rhs; });
  }

  void clear() { data_.clear(); }

  void reserve(size_t size) { data_.reserve(size); }

  void erase(typename std::vector<std::pair<K, V>>::iterator it) { data_.erase(it); }

  auto insert(typename std::vector<std::pair<K, V>>::iterator it,
              std::pair<K, V> value) {
    return data_.insert(it, value);
  }

  void Union(const ContiguousMap& other) {
    std::vector<std::pair<K, V>> result;
    std::set_union(data_.begin(), data_.end(), other.data_.begin(), other.data_.end(),
                   std::back_inserter(result));
    data_ = std::move(result);
  }

  auto Insert(K key, V value) {
    auto it = LowerBound(key);
    if (it != end() and it->first == key) {
      it->second = value;
      return it;
    } else {
      return insert(it, {key, value});
    }
  }

  V& operator[](K key) {
    auto it = LowerBound(key);
    if (it != end() and it->first == key) {
      return it->second;
    } else {
      return insert(it, {key, V{}})->second;
    }
  }

 private:
  ContiguousMap(const ContiguousMap&) = default;
  std::vector<std::pair<K, V>> data_;
};
#pragma once

#include <vector>
#include <algorithm>

template <typename K, typename V>
class ContiguousMap {
 public:
  using value_type = std::pair<K, V>;
  using storage_type = std::vector<value_type>;
  using iterator = typename storage_type::iterator;
  using const_iterator = typename storage_type::const_iterator;

  ContiguousMap() = default;
  ContiguousMap(ContiguousMap&&) = default;
  ContiguousMap& operator=(ContiguousMap&&) = default;
  ContiguousMap& operator=(const ContiguousMap&) = delete;

  ContiguousMap Copy() const { return ContiguousMap{*this}; }

  const_iterator begin() const { return data_.begin(); }

  const_iterator end() const { return data_.end(); }

  const_iterator find(const K& key) const {
    return std::lower_bound(
        data_.begin(), data_.end(), key,
        [](const value_type& lhs, const K& rhs) { return lhs.first < rhs; });
  }

  bool operator==(const ContiguousMap& other) const { return data_ == other.data_; }

  bool operator!=(const ContiguousMap& other) const { return data_ != other.data_; }

  bool operator<(const ContiguousMap& other) const { return data_ < other.data_; }

  bool empty() const { return data_.empty(); }

  size_t size() const { return data_.size(); }

  iterator begin() { return data_.begin(); }

  iterator end() { return data_.end(); }

  iterator find(const K& key) {
    return std::lower_bound(
        data_.begin(), data_.end(), key,
        [](const value_type& lhs, const K& rhs) { return lhs.first < rhs; });
  }

  void clear() { data_.clear(); }

  void reserve(size_t size) { data_.reserve(size); }

  void erase(iterator it) { data_.erase(it); }

  std::pair<iterator, bool> insert(const value_type& value) {
    auto it = find(value.first);
    if (it != end() and it->first == value.first) {
      return {it, false};
    } else {
      return {Insert(it, value), true};
    }
  }

  std::pair<iterator, bool> insert(value_type&& value) {
    auto it = find(value.first);
    if (it != end() and it->first == value.first) {
      return {it, false};
    } else {
      return {Insert(it, std::forward<value_type>(value)), true};
    }
  }

  iterator insert_or_assign(value_type&& value) {
    auto it = find(value.first);
    if (it != end() and it->first == value.first) {
      it->second = std::forward<V>(value.second);
      return it;
    } else {
      return Insert(it, value);
    }
  }

  iterator insert_or_assign(const value_type& value) {
    auto it = find(value.first);
    if (it != end() and it->first == value.first) {
      it->second = value.second;
      return it;
    } else {
      return Insert(it, value);
    }
  }

  V& operator[](const K& key) {
    auto it = find(key);
    if (it != end() and it->first == key) {
      return it->second;
    } else {
      return Insert(it, {key, V{}})->second;
    }
  }

  const V& at(const K& key) const {
    auto it = find(key);
    Assert(it != end());
    return it->second;
  }

  V& at(const K& key) {
    auto it = find(key);
    Assert(it != end());
    return it->second;
  }

  bool Contains(const K& key) const { return find(key) != end(); }

  void Union(const ContiguousMap& other) {
    storage_type result;
    result.reserve(std::max(data_.size(), other.data_.size()));
    std::set_union(data_.begin(), data_.end(), other.data_.begin(), other.data_.end(),
                   std::back_inserter(result));
    data_ = std::move(result);
  }

  iterator Insert(iterator it, value_type&& value) {
    return data_.insert(it, std::forward<value_type>(value));
  }

  iterator Insert(iterator it, const value_type& value) {
    return data_.insert(it, value);
  }

 private:
  ContiguousMap(const ContiguousMap&) = default;
  storage_type data_;
};
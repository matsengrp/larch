#pragma once

#include <vector>
#include <algorithm>

template <typename T>
class ContiguousSet {
 public:
  using storage_type = std::vector<T>;
  using iterator = typename storage_type::iterator;
  using const_iterator = typename storage_type::const_iterator;

  ContiguousSet() = default;
  ContiguousSet(ContiguousSet&&) = default;
  ContiguousSet& operator=(ContiguousSet&&) = default;
  ContiguousSet& operator=(const ContiguousSet&) = delete;

  ContiguousSet Copy() const { return ContiguousSet{*this}; }

  const_iterator begin() const { return data_.begin(); }

  const_iterator end() const { return data_.end(); }

  const_iterator find(const T& value) const {
    return std::lower_bound(data_.begin(), data_.end(), value);
  }

  bool operator==(const ContiguousSet& other) const { return data_ == other.data_; }

  bool operator!=(const ContiguousSet& other) const { return data_ != other.data_; }

  bool operator<(const ContiguousSet& other) const { return data_ < other.data_; }

  bool empty() const { return data_.empty(); }

  size_t size() const { return data_.size(); }

  iterator begin() { return data_.begin(); }

  iterator end() { return data_.end(); }

  iterator find(const T& value) {
    return std::lower_bound(data_.begin(), data_.end(), value);
  }

  void clear() { data_.clear(); }

  void reserve(size_t size) { data_.reserve(size); }

  void erase(iterator it) { data_.erase(it); }

  std::pair<iterator, bool> insert(const T& value) {
    auto it = find(value);
    if (it != data_.end() and *it == value) {
      return {it, false};
    }
    return {data_.insert(it, value), true};
  }

  std::pair<iterator, bool> insert(T&& value) {
    auto it = find(value);
    if (it != data_.end() and *it == value) {
      return {it, false};
    }
    return {data_.insert(it, std::forward<T>(value)), true};
  }

  void Union(const ContiguousSet& other) {
    storage_type result;
    result.reserve(std::max(data_.size(), other.data_.size()));
    ranges::set_union(data_, other.data_, ranges::back_inserter(result));
    data_ = std::move(result);
  }

 private:
  ContiguousSet(const ContiguousSet&) = default;
  storage_type data_;
};

#pragma once

#include <vector>
#include <algorithm>
#include <initializer_list>

template <typename K, typename V>
class ContiguousMap {
 public:
  using value_type = std::pair<K, V>;
  using storage_type = std::vector<value_type>;
  using iterator = typename storage_type::iterator;
  using const_iterator = typename storage_type::const_iterator;

  ContiguousMap() = default;
  ContiguousMap(ContiguousMap&&) noexcept = default;
  ContiguousMap& operator=(ContiguousMap&&) noexcept = default;
  ContiguousMap& operator=(const ContiguousMap&) = delete;
  ~ContiguousMap() = default;

  ContiguousMap(std::initializer_list<value_type> init) : data_{init} {
    data_ |= ranges::actions::sort | ranges::actions::unique;
  }

  ContiguousMap Copy() const { return ContiguousMap{*this}; }

  const_iterator begin() const { return data_.begin(); }

  const_iterator end() const { return data_.end(); }

  const_iterator find(const K& key) const {
    auto result = Find(key);
    if (result != data_.end() and not(result->first == key)) {
      return data_.end();
    }
    return result;
  }

  bool operator==(const ContiguousMap& other) const { return data_ == other.data_; }

  bool operator!=(const ContiguousMap& other) const { return data_ != other.data_; }

  bool operator<(const ContiguousMap& other) const { return data_ < other.data_; }

  bool empty() const { return data_.empty(); }

  size_t size() const { return data_.size(); }

  iterator begin() { return data_.begin(); }

  iterator end() { return data_.end(); }

  iterator find(const K& key) {
    auto result = Find(key);
    if (result != data_.end() and not(result->first == key)) {
      return data_.end();
    }
    return result;
  }

  void clear() { data_.clear(); }

  void reserve(size_t size) { data_.reserve(size); }

  void erase(iterator it) { data_.erase(it); }

  std::pair<iterator, bool> insert(const value_type& value) {
    auto it = Find(value.first);
    if (it != end() and it->first == value.first) {
      return {it, false};
    } else {
      return {Insert(it, value), true};
    }
  }

  std::pair<iterator, bool> insert(value_type&& value) {
    auto it = Find(value.first);
    if (it != end() and it->first == value.first) {
      return {it, false};
    } else {
      return {Insert(it, std::forward<value_type>(value)), true};
    }
  }

  template <typename Mapped>
  std::pair<iterator, bool> insert_or_assign(const K& key, Mapped&& value) {
    auto it = Find(key);
    if (it != end() and it->first == key) {
      it->second = std::forward<Mapped>(value);
      return {it, false};
    } else {
      return {Insert(it, {key, std::forward<Mapped>(value)}), true};
    }
  }

  V& operator[](const K& key) {
    auto it = Find(key);
    if (it != end() and it->first == key) {
      return it->second;
    } else {
      return Insert(it, value_type{key, V{}})->second;
    }
  }

  V& operator[](K&& key) {
    auto it = Find(key);
    if (it != end() and it->first == key) {
      return it->second;
    } else {
      return Insert(it, value_type{std::forward<K>(key), V{}})->second;
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

  bool Contains(const K& key) const { return find(key) != data_.end(); }

  void Union(const ContiguousMap& other) {
    storage_type result;
    result.reserve(std::max(data_.size(), other.data_.size()));
    std::set_union(data_.begin(), data_.end(), other.data_.begin(), other.data_.end(),
                   std::back_inserter(result));
    data_ = std::move(result);
    AssertOrdered();
  }

 private:
  ContiguousMap(const ContiguousMap&) = default;

  void AssertOrdered() const {
#ifndef NDEBUG
    for (size_t i = 0; i + 1 < data_.size(); ++i) {
      Assert((data_[i].first < data_[i + 1].first) or
             (data_[i].first == data_[i + 1].first));
    }
#endif
  }

  iterator Find(const K& key) {
    AssertOrdered();
    return std::lower_bound(
        data_.begin(), data_.end(), key,
        [](const value_type& lhs, const K& rhs) { return lhs.first < rhs; });
  }

  const_iterator Find(const K& key) const {
    AssertOrdered();
    return std::lower_bound(
        data_.begin(), data_.end(), key,
        [](const value_type& lhs, const K& rhs) { return lhs.first < rhs; });
  }

  iterator Insert(iterator it, value_type&& value) {
    auto result = data_.insert(it, std::forward<value_type>(value));
    AssertOrdered();
    return result;
  }

  iterator Insert(iterator it, const value_type& value) {
    auto result = data_.insert(it, value);
    AssertOrdered();
    return result;
  }

  storage_type data_;
};

template <typename K, typename V>
inline std::ostream& operator<<(std::ostream& os, const ContiguousMap<K, V>& map) {
  os << "{ ";
  for (const auto& [key, value] : map) {
    os << key << ": " << value << ", ";
  }
  os << "}";
  return os;
}

#pragma once

#include <vector>
#include <algorithm>

template <typename T, typename Compare = std::less<T>,
          typename Allocator = std::allocator<T>>
class ContiguousSet {
 public:
  using storage_type = std::vector<T, Allocator>;
  using key_type = T;
  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using key_compare = Compare;
  using value_compare = Compare;
  using allocator_type = Allocator;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = typename std::allocator_traits<Allocator>::pointer;
  using const_pointer = typename std::allocator_traits<Allocator>::const_pointer;
  using iterator = typename storage_type::iterator;
  using const_iterator = typename storage_type::const_iterator;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  ContiguousSet() = default;
  ContiguousSet(ContiguousSet&&) noexcept = default;
  ContiguousSet& operator=(ContiguousSet&&) noexcept = default;
  ContiguousSet& operator=(const ContiguousSet&) = delete;
  ~ContiguousSet() = default;

  template <class InputIt>
  ContiguousSet(InputIt first, InputIt last) : data_{first, last} {
    data_ |= ranges::actions::sort(Compare{}) | ranges::actions::unique(Compare{});
  }

  ContiguousSet Copy() const { return ContiguousSet{*this}; }

  const_iterator begin() const { return data_.begin(); }

  const_iterator end() const { return data_.end(); }

  const_iterator find(const T& value) const {
    return std::lower_bound(data_.begin(), data_.end(), value);
  }

  bool Contains(const T& value) const { return find(value) != end(); }

  bool operator==(const ContiguousSet& other) const { return data_ == other.data_; }

  bool operator!=(const ContiguousSet& other) const { return data_ != other.data_; }

  bool operator<(const ContiguousSet& other) const { return data_ < other.data_; }

  bool empty() const { return data_.empty(); }

  size_t size() const { return data_.size(); }

  iterator begin() { return data_.begin(); }

  iterator end() { return data_.end(); }

  iterator find(const T& value) {
    return std::lower_bound(data_.begin(), data_.end(), value, Compare{});
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

  template <typename InputIt>
  void insert(InputIt first, InputIt last) {
    data_.insert(data_.end(), first, last);
    data_ |= ranges::actions::sort(Compare{}) | ranges::actions::unique(Compare{});
  }

  void Union(const ContiguousSet& other) {
    storage_type result;
    result.reserve(std::max(data_.size(), other.data_.size()));
    std::set_union(data_.begin(), data_.end(), other.data_.begin(), other.data_.end(),
                   std::back_inserter(result), Compare{});
    data_ = std::move(result);
  }

 private:
  ContiguousSet(const ContiguousSet&) = default;
  storage_type data_;
};

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const ContiguousSet<T>& set) {
  os << "{ ";
  for (const auto& element : set) {
    os << element << ", ";
  }
  os << "}";
  return os;
}

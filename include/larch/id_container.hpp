#pragma once

#include <vector>
#include <unordered_map>
#include <type_traits>

#include "larch/common.hpp"
#include "larch/tc_vector.hpp"
#include "larch/contiguous_map.hpp"

enum class IdContinuity { Dense, Sparse };

enum class Ordering { Ordered, Unordered };

template <typename Id, typename T, IdContinuity Cont = IdContinuity::Dense,
          Ordering Ord = Ordering::Ordered>
class IdContainer {
  static constexpr auto storage_type_helper = [] {
    if constexpr (Cont == IdContinuity::Dense) {
      if constexpr (std::is_trivially_copyable_v<T>) {
        return type_identity<TCVector<T>>{};
      } else {
        return type_identity<std::vector<T>>{};
      }
    } else {
      if constexpr (Ord == Ordering::Ordered) {
        return type_identity<ContiguousMap<Id, T>>{};
      } else {
        return type_identity<std::unordered_map<Id, T>>{};
      }
    }
  }();

 public:
  static constexpr IdContinuity continuity = Cont;
  static constexpr Ordering ordering =
      Cont == IdContinuity::Dense ? Ordering::Unordered : Ord;

  using storage_type = typename decltype(storage_type_helper)::type;

  using value_type = std::pair<Id, T>;
  using iterator = typename storage_type::iterator;
  using const_iterator = typename storage_type::const_iterator;

  IdContainer() = default;
  IdContainer(IdContainer&&) noexcept = default;
  IdContainer& operator=(IdContainer&&) noexcept = default;
  IdContainer& operator=(const IdContainer&) = delete;
  ~IdContainer() = default;

  IdContainer(std::initializer_list<value_type> init) : data_{init} {}

  IdContainer Copy() const { return IdContainer{*this}; }

  const_iterator begin() const { return data_.begin(); }

  const_iterator end() const { return data_.end(); }

  const_iterator find(Id key) const { return data_.find(key); }

  bool operator==(const IdContainer& other) const { return data_ == other.data_; }

  bool operator!=(const IdContainer& other) const { return data_ != other.data_; }

  bool operator<(const IdContainer& other) const { return data_ < other.data_; }

  bool empty() const { return data_.empty(); }

  size_t size() const { return data_.size(); }

  iterator begin() { return data_.begin(); }

  iterator end() { return data_.end(); }

  iterator find(Id key) { return data_.find(key); }

  void clear() { data_.clear(); }

  void reserve(size_t size) { data_.reserve(size); }

  void erase(iterator it) { data_.erase(it); }

  std::pair<iterator, bool> insert(const value_type& value) {
    return data_.insert(value);
  }

  std::pair<iterator, bool> insert(value_type&& value) {
    if constexpr (Cont == IdContinuity::Dense) {
      bool inserted = value.first.value >= data_.size();
      T& val = this->operator[](value.first);
      val = std::forward<T>(value.second);
      return std::pair<iterator, bool>{std::addressof(val), inserted};
    } else {
      return data_.insert(std::forward<value_type>(value));
    }
  }

  template <typename Mapped>
  std::pair<iterator, bool> insert_or_assign(Id key, Mapped&& value) {
    return data_.insert_or_assign(key, std::forward<Mapped>(value));
  }

  T& operator[](Id key) {
    if constexpr (Cont == IdContinuity::Dense) {
      if (key.value >= data_.size()) {
        data_.resize(key.value + 1);
      }
      return data_[key.value];
    } else {
      return data_[key];
    }
  }

  const T& at(Id key) const {
    if constexpr (Cont == IdContinuity::Dense) {
      return data_.at(key.value);
    } else {
      return data_.at(key);
    }
  }

  T& at(Id key) {
    if constexpr (Cont == IdContinuity::Dense) {
      return data_.at(key.value);
    } else {
      return data_.at(key);
    }
  }

  bool Contains(Id key) const { return data_.find(key) != data_.end(); }

  void Union(const IdContainer& other) { return data_.Union(other); }

  iterator Insert(iterator it, value_type&& value) {
    return data_.Insert(it, std::forward<value_type>(value));
  }

  iterator Insert(iterator it, const value_type& value) {
    return data_.Insert(it, value);
  }

  void resize(size_t size) {
    static_assert(Cont == IdContinuity::Dense);
    data_.resize(size);
  }

  void push_back(T&& value) {
    static_assert(Cont == IdContinuity::Dense);
    data_.push_back(std::forward<T>(value));
  }

 private:
  IdContainer(const IdContainer&) = default;
  storage_type data_;
};

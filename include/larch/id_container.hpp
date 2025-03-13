#pragma once

#include <vector>
#include <unordered_map>
#include <type_traits>
#include <memory>

#include "larch/common.hpp"
#include "larch/contiguous_map.hpp"

static inline constexpr IdContinuity DefIdCont = IdContinuity::Sparse;

template <typename Id, typename T, IdContinuity Cont = DefIdCont,
          Ordering Ord = Ordering::Ordered>
class IdContainer {
 public:
  static constexpr IdContinuity continuity = Cont;
  static constexpr Ordering ordering = Ordering::Unordered;
  // FIXME continuity == IdContinuity::Dense ? Ordering::Unordered : Ord;

 private:
  static constexpr auto storage_type_helper = [] {
    if constexpr (continuity == IdContinuity::Dense) {
      if constexpr (std::is_trivially_copyable_v<T>) {
        return type_identity<std::vector<T>>{};
      } else {
        return type_identity<std::vector<T>>{};
      }
    } else {
      if constexpr (ordering == Ordering::Ordered) {
        return type_identity<ContiguousMap<Id, T>>{};
      } else {
        return type_identity<std::unordered_map<Id, T>>{};
      }
    }
  }();

 public:
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
    if constexpr (continuity == IdContinuity::Dense) {
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
    if constexpr (continuity == IdContinuity::Dense) {
      if (key.value >= data_.size()) {
        data_.resize(key.value + 1);
      }
      return data_[key.value];
    } else {
      return data_[key];
    }
  }

  const T& at(Id key) const {
    if constexpr (continuity == IdContinuity::Dense) {
      return data_.at(key.value);
    } else {
      return data_[key];
    }
  }

  T& at(Id key) {
    if constexpr (continuity == IdContinuity::Dense) {
      return data_.at(key.value);
    } else {
      return data_[key];
    }
  }

  bool Contains(Id key) const {
    if constexpr (continuity == IdContinuity::Dense) {
      return key.value < data_.size();
    } else {
      return data_.find(key) != data_.end();
    }
  }

  Id GetNextAvailableId() const {
    if constexpr (continuity == IdContinuity::Dense) {
      return {data_.size()};
    } else if constexpr (ordering == Ordering::Ordered) {
      if (data_.empty()) {
        return {0};
      }
      return {data_.back().first.value + 1};
    } else {
      auto keys = data_ | ranges::views::keys | ranges::views::filter([](auto& i) {
                    return i.value != MV_UA_NODE_ID;
                  });
      auto result = ranges::max_element(keys);  // TODO linear search
      if (result == keys.end()) {
        return {0};
      }
      return {result->first.value + 1};
    }
  }

  void Union(const IdContainer& other) { return data_.Union(other); }

  iterator Insert(iterator it, value_type&& value) {
    return data_.Insert(it, std::forward<value_type>(value));
  }

  iterator Insert(iterator it, const value_type& value) {
    return data_.Insert(it, value);
  }

  void resize(size_t size) {
    if constexpr (continuity == IdContinuity::Dense) {
      data_.resize(size);
    } else {
      data_.reserve(size);
    }
  }

  void push_back(T&& value) {
    if constexpr (continuity == IdContinuity::Dense) {
      data_.push_back(std::forward<T>(value));
    } else {
      Fail("Can't pushe_back on sparse IDs");
    }
  }

  T* At(Id key) {
    if constexpr (continuity == IdContinuity::Dense) {
      if (key.value >= data_.size()) {
        return nullptr;
      }
      return std::addressof(data_[key.value]);
    } else {
      auto it = data_.find(key);
      if (it == data_.end()) {
        return std::addressof(data_[key]);
      }
      return std::addressof(it->second);
    }
  }

  const T* At(Id key) const {
    if constexpr (continuity == IdContinuity::Dense) {
      if (key.value >= data_.size()) {
        return nullptr;
      }
      return std::addressof(data_[key.value]);
    } else {
      auto it = data_.find(key);
      if (it == data_.end()) {
        Fail("Out of bounds for sparse IDs");
        return nullptr;
      }
      return std::addressof(it->second);
    }
  }

 private:
  IdContainer(const IdContainer&) = default;
  mutable storage_type data_;
};

template <typename Lhs, typename Rhs>
struct ContainerEquivalent : std::false_type {};

template <typename Lhs, typename Rhs, typename Id, IdContinuity Cont, Ordering Ord>
struct ContainerEquivalent<IdContainer<Id, Lhs, Cont, Ord>,
                           IdContainer<Id, Rhs, Cont, Ord>>
    : FeatureEquivalent<Lhs, Rhs> {};

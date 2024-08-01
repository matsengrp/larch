#pragma once

#include <vector>
#include <unordered_map>
#include <type_traits>
#include <memory>

#include "larch/common.hpp"
#include "larch/contiguous_map.hpp"

enum class IdContinuity { Dense, Sparse };

enum class Ordering { Ordered, Unordered };

template <typename Id, typename T, IdContinuity Cont = IdContinuity::Dense,
          Ordering Ord = Ordering::Ordered>
class IdContainer {
  static constexpr auto storage_type_helper = [] {
    if constexpr (Cont == IdContinuity::Dense) {
      if constexpr (std::is_trivially_copyable_v<T>) {
        return type_identity<std::vector<T>>{};
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

  size_t size() const {
    if constexpr (Cont == IdContinuity::Dense) {
      return data_.size();
    } else {
      return std::max(init_size_, data_.size());
    }
  }

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
      if (key.value < init_size_) {
        static const T empty{};
        return empty;
      } else {
        return data_.at(key);
      }
    }
  }

  T& at(Id key) {
    if constexpr (Cont == IdContinuity::Dense) {
      return data_.at(key.value);
    } else {
      if (key.value < init_size_) {
        return data_[key];
      } else {
        return data_.at(key);
      }
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
    if constexpr (Cont == IdContinuity::Dense) {
      data_.resize(size);
    } else {
      init_size_ = size;
    }
  }

  void push_back(T&& value) {
    if constexpr (Cont == IdContinuity::Dense) {
      data_.push_back(std::forward<T>(value));
    } else {
      data_[Id{init_size_++}] = value;
    }
  }

  T* At(Id key) {
    if constexpr (Cont == IdContinuity::Dense) {
      if (key.value >= data_.size()) {
        return nullptr;
      }
      return std::addressof(data_[key.value]);
    } else {
      auto it = data_.find(key);
      if (it == data_.end()) {
        if (key.value < init_size_) {
          return std::addressof(data_[key.value]);
        } else {
          return nullptr;
        }
      }
      return std::addressof(*it);
    }
  }

  const T* At(Id key) const {
    if constexpr (Cont == IdContinuity::Dense) {
      if (key.value >= data_.size()) {
        return nullptr;
      }
      return std::addressof(data_[key.value]);
    } else {
      auto it = data_.find(key);
      if (it == data_.end()) {
        if (key.value >= init_size_) {
          Fail("Out of bounds");
        }
        return nullptr;
      }
      return std::addressof(*it);
    }
  }

 private:
  IdContainer(const IdContainer&) = default;
  storage_type data_;
  size_t init_size_ = 0;
};

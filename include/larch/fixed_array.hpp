#pragma once

#include "larch/common.hpp"

template <typename T>
class FixedArray {
 public:
  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;
  using iterator = pointer;
  using const_iterator = const_pointer;

  explicit FixedArray(size_t size) : size_{size}, data_{new T[size]} {}

  ~FixedArray() { delete[] data_; }

  FixedArray(FixedArray&& other) : size_{other.size_}, data_{other.data_} {
    other.data_ = nullptr;
  }

  FixedArray& operator=(FixedArray&& other) {
    Assert(this != std::addressof(other));
    delete[] data_;
    size_ = other.size_;
    data_ = other.data_;
    other.data_ = nullptr;
    return *this;
  }

  reference at(size_type pos) {
    Assert(pos < size_);
    return data_[pos];
  }
  const_reference at(size_type pos) const {
    Assert(pos < size_);
    return data_[pos];
  }

  reference operator[](size_type pos) noexcept { return data_[pos]; }
  const_reference operator[](size_type pos) const noexcept { return data_[pos]; }

  iterator begin() noexcept { return std::addressof(data_[0]); }
  const_iterator begin() const noexcept { return std::addressof(data_[0]); }
  const_iterator cbegin() const noexcept { return std::addressof(data_[0]); }
  iterator end() noexcept { return std::addressof(data_[size_]); }
  const_iterator end() const noexcept { return std::addressof(data_[size_]); }
  const_iterator cend() const noexcept { return std::addressof(data_[size_]); }

  bool empty() const noexcept { return size_ == 0; }
  size_type size() const noexcept { return size_; }

 private:
  size_t size_;
  T* data_;
};

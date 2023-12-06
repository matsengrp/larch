// Based on https://gist.github.com/loskutov/8c10e34845c71db5177947a0d15ec55f

#pragma once

#include <algorithm>
#include <type_traits>
#include <stdlib.h>
#include <cmath>
#include <malloc.h>

template <class T>
class TCVector {
  std::size_t capacity_, size_;
  T* data_ = nullptr;

  void ensureCapacity(std::size_t desiredCapacity) {
    static_assert(std::is_trivially_copyable_v<T>);
    if (desiredCapacity <= capacity_) return;
    capacity_ = std::max(capacity_ * 2, desiredCapacity);
    T* newData = static_cast<T*>(reallocarray(data_, capacity_, sizeof(T)));
    data_ = newData ? newData : throw std::bad_alloc();
  }

 public:
  TCVector(const TCVector<T>& other) = delete;
  TCVector<T>& operator=(const TCVector<T>& other) = delete;

  TCVector(TCVector<T>&& other) {
    std::swap(capacity_, other.capacity_);
    std::swap(size_, other.size_);
    std::swap(data_, other.data_);
  }

  TCVector<T>& operator=(TCVector<T>&& other) {
    std::swap(capacity_, other.capacity_);
    std::swap(size_, other.size_);
    std::swap(data_, other.data_);
    return *this;
  }

  using iterator = T*;
  using const_iterator = const T*;

  TCVector(std::size_t n = 0)
      : capacity_{std::max(size_t(1), n)}  // default capacity is 1
        ,
        size_{n},
        data_{static_cast<T*>(malloc(sizeof(T) * capacity_))} {
    if (not data_) throw std::bad_alloc{};
  }

  T* begin() const { return data_; }

  T* end() const { return data_ + size_; }

  void push_back(T elem) {
    ensureCapacity(size_ + 1);
    data_[size_++] = elem;
  }

  void resize(size_t s) {
    if (s > size_) {
      ensureCapacity(s);
    }
    size_ = s;
  }

  void clear() { size_ = 0; }

  T popBack() { return data_[--size_]; }

  T& operator[](size_t ind) const { return data_[ind]; }

  T& at(size_t ind) const {
    return (ind < size_) ? data_[ind] : throw std::out_of_range("TCVector::at");
  }

  size_t size() const { return size_; }
  size_t capacity() const { return capacity_; }

  ~TCVector() noexcept { free(data_); }
};

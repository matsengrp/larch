#pragma once

#include <vector>
#include <mutex>

template <typename T>
class concurrent_vector {
 public:
  using value_type = typename std::vector<T>::value_type;
  using allocator_type = typename std::vector<T>::allocator_type;
  using size_type = typename std::vector<T>::size_type;
  using difference_type = typename std::vector<T>::difference_type;
  using reference = typename std::vector<T>::reference;
  using const_reference = typename std::vector<T>::const_reference;
  using pointer = typename std::vector<T>::pointer;
  using const_pointer = typename std::vector<T>::const_pointer;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
  using reverse_iterator = typename std::vector<T>::reverse_iterator;
  using const_reverse_iterator = typename std::vector<T>::const_reverse_iterator;

  iterator begin() noexcept {
    std::unique_lock lock{mtx_};
    return data_.begin();
  }

  const_iterator begin() const noexcept {
    std::unique_lock lock{mtx_};
    return data_.begin();
  }

  const_iterator cbegin() const noexcept {
    std::unique_lock lock{mtx_};
    return data_.cbegin();
  }

  iterator end() noexcept {
    std::unique_lock lock{mtx_};
    return data_.end();
  }

  const_iterator end() const noexcept {
    std::unique_lock lock{mtx_};
    return data_.end();
  }

  const_iterator cend() const noexcept {
    std::unique_lock lock{mtx_};
    return data_.cend();
  }

  size_type size() const noexcept {
    std::unique_lock lock{mtx_};
    return data_.size();
  }

  bool empty() const noexcept {
    std::unique_lock lock{mtx_};
    return data_.empty();
  }

  reference at(size_type pos) {
    std::unique_lock lock{mtx_};
    return data_.at(pos);
  }

  const_reference at(size_type pos) const {
    std::unique_lock lock{mtx_};
    return data_.at(pos);
  }

  iterator push_back(const T& value) {
    std::unique_lock lock{mtx_};
    data_.push_back(value);
    return data_.end() - 1;
  }

  iterator push_back(T&& value) {
    std::unique_lock lock{mtx_};
    data_.push_back(std::move(value));
    return data_.end() - 1;
  }

  void clear() {
    std::unique_lock lock{mtx_};
    data_.clear();
  }

 private:
  std::vector<T> data_;
  mutable std::mutex mtx_;
};

#pragma once

#include <deque>
#include <mutex>

template <typename T>
class concurrent_vector {
 private:
  std::deque<T> data_;
  mutable std::recursive_mutex mtx_;

 public:
   
  auto begin() const noexcept {
    std::lock_guard lock{mtx_};
    return data_.begin();
  }

  auto end() const noexcept {
    std::lock_guard lock{mtx_};
    return data_.end();
  }

  auto begin() noexcept {
    std::lock_guard lock{mtx_};
    return data_.begin();
  }

  auto end() noexcept {
    std::lock_guard lock{mtx_};
    return data_.end();
  }

  size_t size() const noexcept {
    std::lock_guard lock{mtx_};
    return data_.size();
  }

  bool empty() const noexcept {
    std::lock_guard lock{mtx_};
    return data_.empty();
  }

  auto& at(size_t pos) {
    std::lock_guard lock{mtx_};
    return data_.at(pos);
  }

  auto& at(size_t pos) const {
    std::lock_guard lock{mtx_};
    return data_.at(pos);
  }

  decltype(auto) push_back(const T& value) {
    std::lock_guard lock{mtx_};
    data_.push_back(value);
    return data_.end() - 1;
  }

  decltype(auto) push_back(T&& value) {
    std::lock_guard lock{mtx_};
    data_.push_back(std::move(value));
    return data_.end() - 1;
  }

  void clear() {
    std::lock_guard lock{mtx_};
    data_.clear();
  }

};

#pragma once

#include <unordered_map>

#include "larch/parallel/parallel_common.hpp"

template <typename Container>
class Reduction {
  using DataType = std::unordered_map<std::thread::id, Container>;

 public:
  template <typename F, typename... Args>
  void Add(F&& func, Args&&... args) {
    auto rlock = ReadLock(mutex_);
    auto it = data_.find(std::this_thread::get_id());
    if (it == data_.end()) {
      rlock.unlock();
      auto wlock = WriteLock(mutex_);
      Container& container = data_[std::this_thread::get_id()];
      wlock.unlock();
      std::invoke(std::forward<F>(func), static_cast<Container&>(container),
                  std::forward<Args>(args)...);
    } else {
      Container& container = it->second;
      rlock.unlock();
      std::invoke(std::forward<F>(func), static_cast<Container&>(container),
                  std::forward<Args>(args)...);
    }
  }

  template <typename F, typename... Args>
  void Gather(F&& func, Args&&... args) const {
    auto lock = WriteLock(mutex_);
    std::invoke(std::forward<F>(func), static_cast<const DataType&>(data_),
                std::forward<Args>(args)...);
  }

  void clear() {
    auto lock = WriteLock(mutex_);
    data_.clear();
  }

 private:
  mutable std::shared_mutex mutex_;
  DataType data_;
};

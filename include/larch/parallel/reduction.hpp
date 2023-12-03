#pragma once

#include "larch/parallel/parallel_common.hpp"

template <typename Container>
class Reduction {
 public:
  explicit Reduction(size_t buckets) : buckets_{buckets} {}

  template <typename F, typename... Args>
  void AddElement(F&& func, Args&&... args) {
    auto& [mutex, data] = [this]() -> auto& {
      while (true) {
        for (auto& i : buckets_) {
          if (i.mutex_.try_lock()) {
            return i;
          }
        }
      }
    }
    ();
    std::invoke(std::forward<F>(func), static_cast<Container&>(data),
                std::forward<Args>(args)...);
    mutex.unlock();
  }

  template <typename F, typename... Args>
  void GatherAndClear(F&& func, Args&&... args) {
    for (auto& i : buckets_) {
      i.mutex_.lock();
    }
    auto data = buckets_ | ranges::views::transform(
                               [](const auto& i) -> const auto& { return i.data_; });
    std::invoke(std::forward<F>(func), data, std::forward<Args>(args)...);
    buckets_ = FixedArray<Bucket>{buckets_.size()};
    for (auto& i : buckets_) {
      i.mutex_.unlock();
    }
  }

 private:
  struct Bucket {
    std::mutex mutex_;
    Container data_;
  };

  FixedArray<Bucket> buckets_;
};

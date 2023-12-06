#pragma once

#include "larch/parallel/parallel_common.hpp"

template <typename Container>
class Reduction {
 public:
  explicit Reduction(size_t buckets) : buckets_{buckets} {}

  template <typename F, typename... Args>
  decltype(auto) AddElement(F&& func, Args&&... args) {
    auto& [mutex, data] = [this]() -> auto& {
      while (true) {  // TODO endless
        for (auto& i : buckets_) {
          if (i.mutex_.try_lock()) {
            return i;
          }
        }
      }
    }
    ();
    const size_t size = data.size();
    finally cleanup{[this, &mutex, &data, size] {
      size_approx_.fetch_add(data.size() - size);
      mutex.unlock();
    }};
    if constexpr (std::is_void_v<decltype(std::invoke(std::forward<F>(func),
                                                      static_cast<Container&>(data),
                                                      std::forward<Args>(args)...))>) {
      std::invoke(std::forward<F>(func), static_cast<Container&>(data),
                  std::forward<Args>(args)...);
    } else {
      return std::invoke(std::forward<F>(func), static_cast<Container&>(data),
                         std::forward<Args>(args)...);
    }
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
    size_approx_.store(0);
    for (auto& i : buckets_) {
      i.mutex_.unlock();
    }
  }

  size_t size_approx() const { return size_approx_.load(); }

 private:
  struct Bucket {
    std::mutex mutex_;
    Container data_;
  };

  FixedArray<Bucket> buckets_;
  std::atomic<size_t> size_approx_{0};
};

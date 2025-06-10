#pragma once

#include "larch/parallel/parallel_common.hpp"

template <typename Container>
class Reduction {
 public:
  explicit Reduction(size_t buckets)
      : gather_mutex_{}, buckets_{[](size_t num_buckets, std::shared_mutex& mtx) {
          auto lock = WriteLock(mtx);
          return FixedArray<Bucket>{num_buckets};
        }(buckets, gather_mutex_)} {}

  template <typename F, typename... Args>
  decltype(auto) AddElement(F&& func, Args&&... args) {
    auto gather_lock = ReadLock(gather_mutex_);
    auto& [mutex, data] = [this]() -> auto& {
      while (true) {  // TODO endless
        for (auto& i : buckets_) {
          if (i.mutex_.try_lock()) {
            return i;
          }
        }
        std::this_thread::yield();
      }
    }();
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
  decltype(auto) GatherAndClear(F&& func, Args&&... args) {
    auto gather_lock = WriteLock(gather_mutex_);
    finally cleanup{[this] {
      size_approx_.store(0);
      for (auto& i : buckets_) {
        i.data_ = Container{};
      }
    }};
    auto data = buckets_ | ranges::views::transform([](auto& i) -> Container&& {
                  return std::move(i.data_);
                });
    if constexpr (std::is_void_v<decltype(std::invoke(std::forward<F>(func), data,
                                                      std::forward<Args>(args)...))>) {
      std::invoke(std::forward<F>(func), data, std::forward<Args>(args)...);
    } else {
      return std::invoke(std::forward<F>(func), data, std::forward<Args>(args)...);
    }
  }

  size_t size_approx() const { return size_approx_.load(); }

 private:
  struct Bucket {
    std::mutex mutex_;
    Container data_;
  };

  std::shared_mutex gather_mutex_;
  FixedArray<Bucket> buckets_;
  std::atomic<size_t> size_approx_{0};
};

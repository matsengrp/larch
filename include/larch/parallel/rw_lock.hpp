#pragma once

#include "larch/parallel/parallel_common.hpp"

struct rw_spin_impl {
  static constexpr ssize_t READER = 1;
  static constexpr ssize_t WRITER = -1000000;

  /* Call to obtain a reader lock. Must be unlocked later
     with done_read(). Can be upgraded to a writer. */
  inline bool try_read() noexcept {
    if (state_.load() < 0) {
      return false;
    }
    if (state_.fetch_add(+READER) < 0) {
      state_.fetch_add(-READER);
      return false;
    }
    return true;
  }

  /* Call to obtain a writer lock. Must be unlocked later
     with done_write(). Can be downgraded to a reader. */
  inline bool try_write() noexcept {
    ssize_t expected = 0;
    if (not state_.compare_exchange_strong(expected, +WRITER)) {
      return false;
    }
    return true;
  }

  /* Call to upgrade a reader to writer. */
  inline bool try_promote() noexcept {
    ssize_t expected = READER;
    if (not state_.compare_exchange_strong(expected, +WRITER)) {
      return false;
    }
    return true;
  }

  /* Call to downgrade a writer to reader. */
  inline void demote() noexcept { state_.fetch_add(-WRITER + READER); }

  /* Call to unlock a reader. */
  inline void done_read() noexcept { state_.fetch_add(-READER); }

  /* Call to unlock a writer. */
  inline void done_write() noexcept { state_.fetch_add(-WRITER); }

  std::atomic<ssize_t> state_ = 0;
};

class rw_spin_mutex {
 public:
  inline void lock() noexcept {
    while (not impl_.try_write()) {
    }
  }

  inline bool try_lock() noexcept { return impl_.try_write(); }

  inline void unlock() noexcept { impl_.done_write(); }

  inline void lock_shared() noexcept {
    while (not impl_.try_read()) {
    }
  }

  inline bool try_lock_shared() noexcept { return impl_.try_read(); }

  inline void unlock_shared() noexcept { impl_.done_read(); }

 private:
  rw_spin_impl impl_;
};

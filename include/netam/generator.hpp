#pragma once

#include <concepts>
#include <coroutine>
#include <optional>
#include <exception>
#include <iterator>
#include <utility>

template <std::movable T>
class generator {
 public:
  struct promise_type {
    generator<T> get_return_object() {
      return generator{
          std::coroutine_handle<promise_type>::from_promise(*this)};
    }

    static std::suspend_always initial_suspend() noexcept { return {}; }

    static std::suspend_always final_suspend() noexcept { return {}; }

    std::suspend_always yield_value(T value) noexcept {
      value_ = std::move(value);
      return {};
    }

    void await_transform() = delete;

    void unhandled_exception() { exception_ = std::current_exception(); }

    std::optional<T> value_;
    std::exception_ptr exception_;
  };

  using handle_type = std::coroutine_handle<promise_type>;

  generator() = default;
  generator(const generator&) = delete;
  generator& operator=(const generator&) = delete;

  explicit generator(handle_type coroutine) : handle_{coroutine} {}

  generator(generator&& other) noexcept
      : handle_{std::exchange(other.handle_, handle_type{})} {}

  ~generator() {
    if (handle_) {
      handle_.destroy();
    }
  }

  generator& operator=(generator&& other) noexcept {
    if (this != &other) {
      if (handle_) {
        handle_.destroy();
      }
      handle_ = std::exchange(other.handle_, handle_type{});
    }
    return *this;
  }

  class iterator {
   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = const T*;
    using reference = const T&;

    iterator& operator++() {
      handle_.resume();
      return *this;
    }

    void operator++(int) { ++*this; }

    const T& operator*() const {
      const promise_type& promise = handle_.promise();
      if (promise.exception_) {
        std::rethrow_exception(promise.exception_);
      }
      return *promise.value_;
    }

    bool operator==(std::default_sentinel_t) const {
      return (not handle_) or handle_.done();
    }

    explicit iterator(const handle_type coroutine) : handle_{coroutine} {}

   private:
    handle_type handle_;
  };

  iterator begin() {
    if (handle_) {
      handle_.resume();
    }
    return iterator{handle_};
  }

  std::default_sentinel_t end() { return {}; }

 private:
  handle_type handle_;
};

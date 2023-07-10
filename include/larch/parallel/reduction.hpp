#pragma once

#include <vector>
#include <atomic>
#include <shared_mutex>

#include "larch/parallel/node_hashmap.hpp"

template <typename V>
class SizedView {
 public:
  SizedView(V&& view, size_t size) : view_{view}, size_{size} {}
  decltype(auto) begin() { return view_.begin(); }
  decltype(auto) end() { return view_.end(); }
  size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }

 private:
  V view_;
  size_t size_;
};

template <typename T, typename WorkerId = size_t>
class Reduction {
  static constexpr bool UseVector = std::is_same_v<WorkerId, size_t>;
  using Container = std::conditional_t<UseVector, std::vector<std::deque<T>>,
                                       ConcurrentUnorderedMap<WorkerId, std::deque<T>>>;

 public:
  explicit Reduction(size_t workers_count);

  Reduction();
  ~Reduction();

  template <typename... Args>
  std::pair<T&, size_t> Emplace(WorkerId worker, Args&&... args);

  template <typename Lambda>
  void Consume(Lambda&& lambda);

  template <typename Lambda>
  void ConsumeBatches(Lambda&& lambda);

  size_t SizeApprox() const;
  size_t WorkersCount() const;

 private:
  struct Data {
    Container container;
    std::atomic<size_t> size = 0;
  };
  static auto GetRange(Data* data);
  std::atomic<Data*> data_ = nullptr;
  std::shared_mutex mtx_;
};

#include "larch/impl/parallel/reduction_impl.hpp"

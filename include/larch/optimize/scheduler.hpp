#pragma once

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <memory>
#include <atomic>
#include <functional>
#include <set>
#include <map>
#include <type_traits>

class TaskBase {
 public:
  virtual ~TaskBase() {}

 protected:
  friend class Scheduler;
  virtual void Run(size_t worker) = 0;
  virtual bool CanIterate() const = 0;
  virtual bool Finish() = 0;
};

inline bool operator<(const std::reference_wrapper<TaskBase>& lhs,
                      const std::reference_wrapper<TaskBase>& rhs) {
  return std::addressof(lhs.get()) < std::addressof(rhs.get());
}

template <typename F>
class Task : public TaskBase {
 public:
  Task(F&& func);

  void Join();

 private:
  void Run(size_t worker) override;
  bool CanIterate() const override;
  bool Finish() override;

  F func_;
  std::atomic<size_t> iteration_;
  std::atomic<bool> can_iterate_;
  std::mutex mtx_;
  std::condition_variable finished_;
  bool is_finished_ = false;
};

class Scheduler {
 public:
  inline explicit Scheduler(size_t workers_count = std::thread::hardware_concurrency());
  inline ~Scheduler();

  Scheduler(const Scheduler&) = delete;
  Scheduler(Scheduler&&) = delete;
  Scheduler& operator=(const Scheduler&) = delete;
  Scheduler& operator=(Scheduler&&) = delete;

  inline void AddTask(TaskBase& task);

  inline size_t WorkersCount() const;

 private:
  static inline void Worker(Scheduler& self, size_t id);
  const size_t workers_count_;
  std::vector<std::thread> workers_;
  std::atomic<bool> destroy_ = false;
  std::mutex mtx_;
  std::deque<std::reference_wrapper<TaskBase>> queue_;
  std::set<std::reference_wrapper<TaskBase>> finished_tasks_;  // TODO: slow
  std::condition_variable queue_not_empty_;
};

template <typename V>
class SizedView {
 public:
  SizedView(V&& view, size_t size) : view_{view}, size_{size} {}
  decltype(auto) begin() { return view_.begin(); }
  decltype(auto) end() { return view_.end(); }
  size_t size() const { return size_; }

 private:
  V view_;
  size_t size_;
};

template <typename T, typename WorkerId = size_t>
class Reduction {
  static constexpr bool UseVector = std::is_same_v<WorkerId, size_t>;
  using Container = std::conditional_t<UseVector, std::vector<std::deque<T>>,
                                       std::map<WorkerId, std::deque<T>>>;

 public:
  explicit Reduction(size_t workers_count) : size_{0} {
    if constexpr (UseVector) {
      data_.resize(workers_count);
    }
  }

  Reduction() : size_{0} { static_assert(not UseVector); }

  template <typename... Args>
  T& Emplace(WorkerId worker, Args&&... args) {
    size_.fetch_add(1);
    return data_[worker].emplace_back(std::forward<Args>(args)...);
  }

  auto Get() {
    if constexpr (UseVector) {
      return data_ | ranges::views::all;
    } else {
      return data_ | ranges::views::values;
    }
  }

  auto GetAll() { return SizedView{Get() | ranges::views::join, Size()}; }

  size_t Size() const { return size_.load(); }
  bool Empty() const { return size_.load() > 0; }
  size_t WorkersCount() const { return data_.size(); }

  void Clear() {
    size_.store(0);
    data_.clear();
  }

 private:
  Container data_;
  std::atomic<size_t> size_;
};

template <typename T>
class Snapshot {
 public:
  template <typename... Args>
  explicit Snapshot(Args&&... args) : data_{new T{std::forward<Args>(args)...}} {}

  ~Snapshot() { delete *data_.load(); }

  T& Get() { return *data_.load(); }

  template <typename Lambda, typename... Args>
  void Take(Lambda&& lambda, Args&&... args) {
    T* next = new T{std::forward<Args>(args)...};
    std::unique_ptr<T> curr{std::atomic_exchange(data_, next)};
    lambda(*curr);
  }

 private:
  std::atomic<T*> data_;
};

#include "larch/impl/optimize/scheduler_impl.hpp"

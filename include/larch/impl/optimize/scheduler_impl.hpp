
template <typename F>
Task<F>::Task(F&& func)
    : func_{std::forward<F>(func)}, iteration_{0}, can_iterate_{false} {}

template <typename F>
void Task<F>::Join() {
  std::unique_lock lock{mtx_};
  while (not is_finished_) {
    finished_.wait(lock);
  }
}

template <typename F>
void Task<F>::Run(size_t worker) {
  if (not func_(iteration_.fetch_add(1), worker)) {
    can_iterate_.store(false);
  }
}

template <typename F>
bool Task<F>::CanIterate() const {
  return can_iterate_.load();
}

template <typename F>
bool Task<F>::Finish() {
  std::unique_lock lock{mtx_};
  if (not is_finished_) {
    is_finished_ = true;
    finished_.notify_all();
    return true;
  }
  return false;
}

Scheduler::Scheduler(size_t workers_count) : workers_count_{workers_count} {
  for (size_t i = 0; i < workers_count; ++i) {
    workers_.push_back(std::thread(Worker, std::ref(*this), size_t{i}));
  }
}

Scheduler::~Scheduler() {
  destroy_.store(true);
  {
    std::unique_lock lock{mtx_};
    queue_.clear();
    queue_not_empty_.notify_all();
  }
  for (auto& i : workers_) {
    i.join();
  }
}

void Scheduler::AddTask(TaskBase& task) {
  std::unique_lock lock{mtx_};
  queue_.push_back(std::ref(task));
  queue_not_empty_.notify_all();
}

size_t Scheduler::WorkersCount() const { return workers_count_; }

void Scheduler::Worker(Scheduler& self, size_t id) {
  while (not self.destroy_.load()) {
    std::unique_lock lock{self.mtx_};
    for (auto i : self.finished_tasks_) {
      i.get().Finish();
    }
    self.finished_tasks_.clear();
    while (self.queue_.empty() and not self.destroy_.load()) {
      self.queue_not_empty_.wait(lock);
    }
    if (self.destroy_.load()) {
      return;
    }
    auto& task = self.queue_.front().get();
    if (not task.CanIterate()) {
      if (self.finished_tasks_.insert(std::ref(task)).second) {
        self.queue_.pop_front();
      }
    } else {
      lock.unlock();
      task.Run(id);
    }
  }
}

template <typename T, typename WorkerId>
Reduction<T, WorkerId>::Reduction(size_t workers_count) : size_{0} {
  if constexpr (UseVector) {
    data_.resize(workers_count);
  }
}

template <typename T, typename WorkerId>
Reduction<T, WorkerId>::Reduction() : size_{0} {
  static_assert(not UseVector);
}

template <typename T, typename WorkerId>
template <typename... Args>
T& Reduction<T, WorkerId>::Emplace(WorkerId worker, Args&&... args) {
  size_.fetch_add(1);
#ifdef USE_TSAN
  std::unique_lock lock{mtx_};
#endif
  return data_[worker].emplace_back(std::forward<Args>(args)...);
}

template <typename T, typename WorkerId>
auto Reduction<T, WorkerId>::Get() {
  if constexpr (UseVector) {
    return data_ | ranges::views::all;
  } else {
    return data_ | ranges::views::values;
  }
}

template <typename T, typename WorkerId>
auto Reduction<T, WorkerId>::GetAll() {
#ifdef USE_TSAN
  std::unique_lock lock{mtx_};
#endif
  return SizedView{Get() | ranges::views::join, Size()};
}

template <typename T, typename WorkerId>
size_t Reduction<T, WorkerId>::Size() const {
  return size_.load();
}

template <typename T, typename WorkerId>
bool Reduction<T, WorkerId>::Empty() const {
  return size_.load() > 0;
}

template <typename T, typename WorkerId>
size_t Reduction<T, WorkerId>::WorkersCount() const {
  return data_.size();
}

template <typename T, typename WorkerId>
void Reduction<T, WorkerId>::Clear() {
  size_.store(0);
#ifdef USE_TSAN
  std::unique_lock lock{mtx_};
#endif
  data_.clear();
}

template <typename T>
template <typename... Args>
Snapshot<T>::Snapshot(Args&&... args) : data_{new T{std::forward<Args>(args)...}} {}

template <typename T>
Snapshot<T>::~Snapshot() {
  delete data_.load();
}

template <typename T>
T& Snapshot<T>::Get() {
  return *data_.load();
}

template <typename T>
template <typename Lambda, typename... Args>
void Snapshot<T>::Take(Lambda&& lambda, Args&&... args) {
  T* next = new T{std::forward<Args>(args)...};
  std::unique_ptr<T> curr{std::atomic_exchange(&data_, next)};
  lambda(*curr);
}

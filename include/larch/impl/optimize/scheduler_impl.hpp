
template <typename F>
Task<F>::Task(F&& func) : func_{std::forward<F>(func)}, iteration_{0} {}

template <typename F>
void Task<F>::Join() {
  std::unique_lock lock{mtx_};
  while (not is_finished_) {
    finished_.wait(lock);
  }
}

template <typename F>
bool Task<F>::Run(size_t worker) {
  return func_(iteration_.fetch_add(1), worker);
}

template <typename F>
void Task<F>::Finish() {
  if (not is_finished_) {
    std::unique_lock lock{mtx_};
    if (not is_finished_) {
      is_finished_ = true;
      finished_.notify_all();
    }
  }
}

Scheduler::Scheduler(size_t workers_count) : workers_count_{workers_count} {
  for (size_t i = 0; i < workers_count; ++i) {
    workers_.push_back(std::thread(&Scheduler::Worker, this, size_t{i}));
    running_tasks_.push_back(NoId);
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
  queue_.push_back({std::ref(task), ids_.fetch_add(1)});
  queue_not_empty_.notify_all();
}

size_t Scheduler::WorkersCount() const { return workers_count_; }

void Scheduler::Worker(size_t id) {
  while (not destroy_.load()) {
    std::unique_lock lock{mtx_};
    while (queue_.empty() and not destroy_.load()) {
      queue_not_empty_.wait(lock);
    }
    if (destroy_.load()) {
      return;
    }
    running_tasks_.at(id) = NoId;
    bool front_task_done = false;
    while (not done_tasks_.empty() and not queue_.empty()) {
      size_t task_id = queue_.front().second;
      auto task = done_tasks_.find(task_id);
      if (task != done_tasks_.end()) {
        front_task_done = true;
        if (not IsTaskRunning(task_id)) {
          queue_.front().first.get().Finish();
          queue_.pop_front();
          done_tasks_.erase(task);
          front_task_done = false;
          continue;
        }
      }
      break;
    }
    if (not queue_.empty() and not front_task_done) {
      QueueItem task = queue_.front();
      running_tasks_.at(id) = task.second;
      lock.unlock();
      if (not task.first.get().Run(id)) {
        lock.lock();
        running_tasks_.at(id) = NoId;
        done_tasks_.insert(task.second);
      }
    }
  }
}

bool Scheduler::IsTaskRunning(size_t task_id) const {
  return std::find(running_tasks_.begin(), running_tasks_.end(), task_id) !=
         running_tasks_.end();
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

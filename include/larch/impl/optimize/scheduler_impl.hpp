
template <typename F>
Task<F>::Task(F&& func) : func_{std::forward<F>(func)} {}

template <typename F>
void Task<F>::Join() {
  std::unique_lock lock{mtx_};
  while (not is_finished_) {
    finished_.wait(lock);
  }
}

template <typename F>
void Task<F>::Run() {
  if (not func_(iteration_.fetch_add(1))) {
    can_iterate_.store(false);
  }
}

template <typename F>
bool Task<F>::CanIterate() const {
  return can_iterate_.load();
}

template <typename F>
void Task<F>::Finish() {
  std::unique_lock lock{mtx_};
  if (not is_finished_) {
    is_finished_ = true;
    finished_.notify_all();
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

void Scheduler::Start() {
  size_t thread_count = std::thread::hardware_concurrency();
  for (size_t i = 0; i < thread_count; ++i) {
    workers_.push_back(std::thread(Worker, std::ref(*this)));
  }
}

void Scheduler::AddTask(TaskBase& task) {
  std::unique_lock lock{mtx_};
  queue_.push_back(std::ref(task));
  queue_not_empty_.notify_one();
}

void Scheduler::Worker(Scheduler& self) {
  while (not self.destroy_.load()) {
    std::unique_lock lock{self.mtx_};
    while (self.queue_.empty() and not self.destroy_.load()) {
      self.queue_not_empty_.wait(lock);
    }
    if (self.destroy_.load()) {
      return;
    }
    auto& task = self.queue_.front().get();
    if (not task.CanIterate()) {
      self.queue_.pop_front();
      task.Finish();
    } else {
      lock.unlock();
      task.Run();
    }
  }
}

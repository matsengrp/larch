TaskBase::TaskBase(Scheduler& scheduler) : id_{scheduler.NewTaskId()} {}

size_t TaskBase::GetId() const { return id_; }

template <typename F>
Task<F>::Task(Scheduler& scheduler, F&& func)
    : TaskBase{scheduler}, func_{std::forward<F>(func)} {}

template <typename F>
Task<F>::~Task() {
  Join();
}

template <typename F>
void Task<F>::Join() {
  std::unique_lock lock{done_mtx_};
  while (not done_) {
    is_done_.wait(lock);
  }
}

template <typename F>
void Task<F>::Join(Scheduler& scheduler) {
  if (not scheduler.JoinTask(*this)) {
    Join();
  }
}

template <typename F>
bool Task<F>::Run(size_t worker) {
  return func_(iteration_.fetch_add(1), worker);
}

template <typename F>
bool Task<F>::Finish(size_t workers_count) {
  if (finished_.fetch_add(1) + 1 == workers_count) {
    std::unique_lock lock{done_mtx_};
    done_ = true;
    is_done_.notify_all();
    return true;
  }
  return false;
}

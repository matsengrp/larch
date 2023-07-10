
template <typename T, typename WorkerId>
Reduction<T, WorkerId>::Reduction(size_t workers_count) : data_{new Data} {
  static_assert(UseVector);
  data_.load()->container.resize(workers_count);
}

template <typename T, typename WorkerId>
Reduction<T, WorkerId>::Reduction() : data_{new Data} {
  static_assert(not UseVector);
}

template <typename T, typename WorkerId>
Reduction<T, WorkerId>::~Reduction() {
  delete data_.exchange(nullptr);
}

template <typename T, typename WorkerId>
template <typename... Args>
std::pair<T&, size_t> Reduction<T, WorkerId>::Emplace(WorkerId worker, Args&&... args) {
  std::shared_lock lock{mtx_};  // Intentionally reversed role (writer = shared)
  Data* data = data_.load();
  size_t size = data->size.fetch_add(1) + 1;
  if constexpr (UseVector) {
    return {data->container[worker].emplace_back(std::forward<Args>(args)...), size};
  } else {
    return {data->container.AtDefault(worker).GetExclusive(
                [&](auto& val, auto&&... lambda_args) -> decltype(auto) {
                  return val.emplace_back(
                      std::forward<decltype(lambda_args)>(lambda_args)...);
                },
                std::forward<Args>(args)...),
            size};
  }
}

template <typename T, typename WorkerId>
template <typename Lambda>
void Reduction<T, WorkerId>::Consume(Lambda&& lambda) {
  std::unique_lock lock{mtx_};
  Data* old = data_.exchange(new Data);
  lambda(SizedView{GetRange(old) | ranges::views::join, old->size.load()});
  delete old;
}

template <typename T, typename WorkerId>
template <typename Lambda>
void Reduction<T, WorkerId>::ConsumeBatches(Lambda&& lambda) {
  std::unique_lock lock{mtx_};
  Data* old = data_.exchange(new Data);
  lambda(GetRange(old));
  delete old;
}

template <typename T, typename WorkerId>
size_t Reduction<T, WorkerId>::SizeApprox() const {
  Data* data = data_.load();
  return data->size.load();
}

template <typename T, typename WorkerId>
size_t Reduction<T, WorkerId>::WorkersCount() const {
  Data* data = data_.load();
  return data->container.size();
}

template <typename T, typename WorkerId>
auto Reduction<T, WorkerId>::GetRange(Data* data) {
  if constexpr (UseVector) {
    return data->container | ranges::views::all;
  } else {
    return data->container.All() | ranges::views::values;
  }
}

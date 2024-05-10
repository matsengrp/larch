#pragma once

#include <chrono>

class Benchmark {
 public:
  using TimePoint = decltype(std::chrono::high_resolution_clock::now());

  inline Benchmark(bool start_on_init = true);

  inline void start();
  inline void stop();

  template <typename time_scale>
  inline auto lap();
  inline auto lapMs();
  inline std::string lapFormat();

  template <typename time_scale>
  inline auto duration() const;
  inline auto durationUs() const;
  inline auto durationMs() const;
  inline auto durationS() const;
  inline std::string durationFormat() const;

  inline static std::string formatMs(long int ms_count);

 private:
  TimePoint start_;
  TimePoint stop_;
};

///////////////////////////////////////////////////////////////////////////////

Benchmark::Benchmark(bool start_on_init) {
  if (start_on_init) {
    start();
  }
}

void Benchmark::start() { start_ = std::chrono::high_resolution_clock::now(); }

void Benchmark::stop() { stop_ = std::chrono::high_resolution_clock::now(); }

template <typename time_scale>
auto Benchmark::duration() const {
  return std::chrono::duration_cast<time_scale>(stop_ - start_).count();
}

auto Benchmark::durationUs() const { return duration<std::chrono::microseconds>(); }

auto Benchmark::durationMs() const { return duration<std::chrono::milliseconds>(); }

auto Benchmark::durationS() const { return duration<std::chrono::seconds>(); }

std::string Benchmark::durationFormat() const { return formatMs(durationMs()); }

template <typename time_scale>
auto Benchmark::lap() {
  stop();
  auto result = duration<time_scale>();
  start_ = stop_ = std::chrono::high_resolution_clock::now();
  return result;
}

auto Benchmark::lapMs() { return lap<std::chrono::milliseconds>(); }

std::string Benchmark::lapFormat() { return formatMs(lapMs()); }

std::string Benchmark::formatMs(long int ms_count) {
  std::stringstream ss;
  auto min = ms_count / (60 * 1000);
  auto sec = static_cast<double>(ms_count % (60 * 1000)) / 1000.0;
  ss << min << "m" << sec << "s";
  return ss.str();
}

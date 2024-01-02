#pragma once

#include <chrono>

class Benchmark {
 public:
  using TimePoint = decltype(std::chrono::high_resolution_clock::now());

  inline void start();
  inline void stop();

  inline auto lapMs();

  inline auto durationUs() const;
  inline auto durationMs() const;
  inline auto durationS() const;

 private:
  TimePoint start_;
  TimePoint stop_;
};

///////////////////////////////////////////////////////////////////////////////

void Benchmark::start() { start_ = std::chrono::high_resolution_clock::now(); }

void Benchmark::stop() { stop_ = std::chrono::high_resolution_clock::now(); }

auto Benchmark::durationUs() const {
  return std::chrono::duration_cast<std::chrono::microseconds>(stop_ - start_).count();
}

auto Benchmark::durationMs() const {
  return std::chrono::duration_cast<std::chrono::milliseconds>(stop_ - start_).count();
}

auto Benchmark::durationS() const {
  return std::chrono::duration_cast<std::chrono::seconds>(stop_ - start_).count();
}

auto Benchmark::lapMs() {
  stop();
  auto result = durationMs();
  start_ = stop_ = std::chrono::high_resolution_clock::now();
  return result;
}
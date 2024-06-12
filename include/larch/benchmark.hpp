#pragma once

#include <chrono>

// template <typename time_scale>
class Benchmark {
 public:
  using TimePoint = decltype(std::chrono::high_resolution_clock::now());
  using Us = std::chrono::microseconds;
  using Ms = std::chrono::milliseconds;
  using S = std::chrono::seconds;

  inline Benchmark(bool start_on_init = true);

  inline void start();
  inline void stop();

  template <typename time_scale>
  inline auto lap();
  inline auto lapUs();
  inline auto lapMs();
  inline auto lapS();

  template <typename time_scale>
  inline std::string lapFormat();
  inline std::string lapFormatUs();
  inline std::string lapFormatMs();
  inline std::string lapFormatS();

  template <typename time_scale>
  inline auto duration() const;
  inline auto durationUs() const;
  inline auto durationMs() const;
  inline auto durationS() const;

  template <typename time_scale>
  inline std::string durationFormat() const;
  inline std::string durationFormatUs() const;
  inline std::string durationFormatMs() const;
  inline std::string durationFormatS() const;

  template <typename time_scale>
  inline static std::string format(long int ticks);
  inline static std::string formatUs(long int ticks);
  inline static std::string formatMs(long int ticks);
  inline static std::string formatS(long int ticks);

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

template <typename time_scale>
std::string Benchmark::durationFormat() const {
  return format<time_scale>(duration<time_scale>());
}

std::string Benchmark::durationFormatUs() const {
  return durationFormat<std::chrono::microseconds>();
}

std::string Benchmark::durationFormatMs() const {
  return durationFormat<std::chrono::milliseconds>();
}

std::string Benchmark::durationFormatS() const {
  return durationFormat<std::chrono::seconds>();
}

template <typename time_scale>
auto Benchmark::lap() {
  stop();
  auto result = duration<time_scale>();
  start_ = stop_ = std::chrono::high_resolution_clock::now();
  return result;
}

auto Benchmark::lapUs() { return lap<std::chrono::microseconds>(); }

auto Benchmark::lapMs() { return lap<std::chrono::milliseconds>(); }

auto Benchmark::lapS() { return lap<std::chrono::seconds>(); }

template <typename time_scale>
std::string Benchmark::lapFormat() {
  return format<time_scale>(lap<time_scale>());
}

std::string Benchmark::lapFormatUs() { return lapFormat<std::chrono::microseconds>(); }

std::string Benchmark::lapFormatMs() { return lapFormat<std::chrono::milliseconds>(); }

std::string Benchmark::lapFormatS() { return lapFormat<std::chrono::seconds>(); }

template <typename time_scale>
std::string Benchmark::format(long int ticks) {
  long int ticks_per_second;
  if constexpr (std::is_same_v<time_scale, std::chrono::seconds>) {
    ticks_per_second = 1;
  } else if constexpr (std::is_same_v<time_scale, std::chrono::milliseconds>) {
    ticks_per_second = 1000;
  } else if constexpr (std::is_same_v<time_scale, std::chrono::microseconds>) {
    ticks_per_second = 1000000;
  } else {
    static_assert(!std::is_same_v<time_scale, time_scale>, "ERROR: Unsupported type.");
  }

  std::stringstream ss;
  auto min = ticks / (60 * ticks_per_second);
  auto sec = static_cast<double>(ticks % (60 * ticks_per_second)) /
             static_cast<double>(ticks_per_second);
  ss << min << "m" << sec << "s";
  return ss.str();
}

std::string Benchmark::formatUs(long int ticks) {
  return format<std::chrono::microseconds>(ticks);
}

std::string Benchmark::formatMs(long int ticks) {
  return format<std::chrono::milliseconds>(ticks);
}

std::string Benchmark::formatS(long int ticks) {
  return format<std::chrono::seconds>(ticks);
}

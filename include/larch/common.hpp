#pragma once

#include <cstddef>
#include <limits>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wstack-usage="
#include <range/v3/action/push_back.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/action/unique.hpp>
#include <range/v3/algorithm/unique.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/counted.hpp>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/group_by.hpp>
#include <range/v3/view/indices.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/subrange.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>
#pragma GCC diagnostic pop

struct NodeId;
struct EdgeId;
struct CladeIdx;

static constexpr const size_t NoId = std::numeric_limits<size_t>::max();

template <typename T, typename Id>
[[nodiscard]] static T& GetOrInsert(std::vector<T>& data, Id id) {
  size_t idx{};
  if constexpr (std::is_same_v<Id, size_t>) {
    idx = id;
  } else {
    idx = id.value;
  }

  if (idx >= data.size()) {
    data.resize(idx + 1);
  }
  return data[idx];
}

template <typename T>
void MoveElements(std::vector<T>&& source, std::vector<T>& destination) {
  destination.resize(0);
  destination.insert(destination.end(), std::make_move_iterator(source.begin()),
                     std::make_move_iterator(source.end()));
  source.resize(0);
}

namespace Transform {
template <typename T>
inline auto To() {
  return ranges::views::transform([&](auto&& i) { return T{i}; });
}
}  // namespace Transform

inline constexpr const auto HashCombine = [](size_t lhs, size_t rhs) noexcept {
  constexpr const size_t GoldenRatioFractional = 0x9e3779b97f4a7c15;
  constexpr const size_t LeftShift64 = 12;
  constexpr const size_t RightShift64 = 4;
  lhs ^= rhs + GoldenRatioFractional + (lhs << LeftShift64) + (lhs >> RightShift64);
  return lhs;
};

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifndef NDEBUG
#include <csignal>
#include <iostream>
#define Assert(x)                                           \
  {                                                         \
    if (not(x)) {                                           \
      std::cerr << "Assert failed: \"" #x "\" in " __FILE__ \
                   ":" TOSTRING(__LINE__) "\n";             \
      std::raise(SIGTRAP);                                  \
    }                                                       \
  }
#else
#define Assert(x)                                                       \
  {                                                                     \
    if (not(x))                                                         \
      throw std::runtime_error("Assert failed: \"" #x "\" in " __FILE__ \
                               ":" TOSTRING(__LINE__));                 \
  }
#endif

[[noreturn]] inline void Fail(const char* msg) { throw std::runtime_error(msg); }

#define MOVE_ONLY(x)                    \
  x(x&&) noexcept = default;            \
  x(const x&) = delete;                 \
  x& operator=(x&&) noexcept = default; \
  x& operator=(const x&) = delete;      \
  ~x() = default

//////////////////////////////////////////////////////////////////////////////////////

template <typename, template <typename...> typename>
struct is_specialization : std::false_type {};

template <template <typename...> typename Template, typename... Args>
struct is_specialization<Template<Args...>, Template> : std::true_type {};

template <typename L, template <typename...> typename R>
constexpr bool is_specialization_v = is_specialization<L, R>::value;

//////////////////////////////////////////////////////////////////////////////////////

template <typename, typename>
struct tuple_contains {};

template <typename... Types, typename Type>
struct tuple_contains<std::tuple<Types...>, Type>
    : std::bool_constant<(std::is_same_v<Types, Type> || ...)> {};

template <typename Tuple, typename Type>
inline constexpr bool tuple_contains_v = tuple_contains<Tuple, Type>::value;

//////////////////////////////////////////////////////////////////////////////////////
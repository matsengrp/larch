#pragma once

#include <cstddef>
#include <limits>
#include <vector>
#include <set>
#include <unordered_set>
#include <tuple>
#include <typeinfo>

//////////////////////////////////////////////////////////////////////////////////////

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wstack-usage="
#include <range/v3/all.hpp>
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
#define Assert(x)                                                       \
  {                                                                     \
    if (not(x)) {                                                       \
      std::cerr << "Assert failed: \"" #x "\" in " __FILE__             \
                   ":" TOSTRING(__LINE__) "\n";                         \
      /*std::raise(SIGTRAP);*/                                          \
      throw std::runtime_error("Assert failed: \"" #x "\" in " __FILE__ \
                               ":" TOSTRING(__LINE__));                 \
    }                                                                   \
  }
#else
#define Assert(x)                                                       \
  {                                                                     \
    if (not(x))                                                         \
      throw std::runtime_error("Assert failed: \"" #x "\" in " __FILE__ \
                               ":" TOSTRING(__LINE__));                 \
  }
#endif

[[noreturn]] inline void Fail(const char* msg) {
#ifndef NDEBUG
  std::cerr << msg << "\n";
  std::raise(SIGTRAP);
#endif
  throw std::runtime_error(msg);
}

#define MOVE_ONLY(x)                    \
  x(x&&) noexcept = default;            \
  x(const x&) = delete;                 \
  x& operator=(x&&) noexcept = default; \
  x& operator=(const x&) = delete;      \
  ~x() = default

#define MOVE_ONLY_VIRT_DTOR(x)          \
  x(x&&) noexcept = default;            \
  x(const x&) = delete;                 \
  x& operator=(x&&) noexcept = default; \
  x& operator=(const x&) = delete;      \
  virtual ~x() = default

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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wc++20-extensions"
template <typename T>
inline const auto tuple_to_string_impl = [
]<std::size_t... I>(std::index_sequence<I...>) {
  std::string result = "std::tuple<";
  result += (... + (std::string{typeid(std::tuple_element_t<I, T>).name()} +
                    std::string{", "}));
  result.pop_back();
  result.pop_back();
  result += ">";
  return result;
};

template <typename T>
inline std::string tuple_to_string() {
  return tuple_to_string_impl<T>(std::make_index_sequence<std::tuple_size_v<T>>());
}
#pragma GCC diagnostic pop

//////////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& vector) {
  os << "[ ";
  for (const auto& element : vector) {
    os << element << " ";
  }
  os << "]";
  return os;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::set<T>& set) {
  os << "{ ";
  for (const auto& element : set) {
    os << element << " ";
  }
  os << "}";
  return os;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::unordered_set<T>& set) {
  os << "{ ";
  for (const auto& element : set) {
    os << element << " ";
  }
  os << "}";
  return os;
}

template <typename K, typename V>
inline std::ostream& operator<<(std::ostream& os, const std::map<K, V>& map) {
  os << "{ ";
  for (const auto& [key, value] : map) {
    os << "( " << key << ": " << value << " ) ";
  }
  os << "}";
  return os;
}

template <typename K, typename V>
inline std::ostream& operator<<(std::ostream& os, const std::unordered_map<K, V>& map) {
  os << "{ ";
  for (const auto& [key, value] : map) {
    os << "( " << key << ": " << value << " ) ";
  }
  os << "}";
  return os;
}

template <typename T1, typename T2>
inline std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& tup) {
  os << "[ " << tup.first << " " << tup.second << "]";
  return os;
}

template <typename T1, typename T2, typename T3>
inline std::ostream& operator<<(std::ostream& os, const std::tuple<T1, T2, T3>& tup) {
  os << "[ " << std::get<0>(tup) << " " << std::get<1>(tup) << " " << std::get<2>(tup)
     << "]";
  return os;
}

//////////////////////////////////////////////////////////////////////////////////////

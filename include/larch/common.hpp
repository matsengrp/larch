#pragma once

#include <cstddef>
#include <limits>
#include <vector>
#include <type_traits>
#include <tuple>

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
inline auto GetParent() {
  return ranges::views::transform([](auto&& i) { return i.GetParent(); });
}
inline auto GetChild() {
  return ranges::views::transform([](auto&& i) { return i.GetChild(); });
}
inline auto GetId() {
  return ranges::views::transform([](auto&& i) { return i.GetId(); });
}
template <typename DAG>
inline auto ToNodes(DAG dag) {
  return ranges::views::transform([dag](auto&& i) {
    return typename DAG::Node{dag, i};
  });
}
template <typename DAG>
inline auto ToEdges(DAG dag) {
  return ranges::views::transform([dag](auto&& i) {
    return typename DAG::Edge{dag, i};
  });
}

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

#define Assert(x)                                                       \
  {                                                                     \
    if (not(x))                                                         \
      throw std::runtime_error("Assert failed: \"" #x "\" in " __FILE__ \
                               ":" TOSTRING(__LINE__));                 \
  }

[[noreturn]] inline void Fail(const char* msg) { throw std::runtime_error(msg); }

#define MOVE_ONLY(x)                    \
  x(x&&) noexcept = default;            \
  x(const x&) = delete;                 \
  x& operator=(x&&) noexcept = default; \
  x& operator=(const x&) = delete;      \
  ~x() = default

// tuple_contians

template <typename, typename>
struct tuple_contians {};

template <typename... Types, typename Type>
struct tuple_contians<std::tuple<Types...>, Type>
    : std::bool_constant<(std::is_same_v<Types, Type> || ...)> {};

template <typename Tuple, typename Type>
inline constexpr bool tuple_contians_v = tuple_contians<Tuple, Type>::value;

// filter

template <template <typename> typename Predicate,
          template <typename...> typename Result, typename... Ts>
struct filter;

template <template <typename> typename Predicate,
          template <typename...> typename Result, typename... Ts>
using filter_t = typename filter<Predicate, Result, Ts...>::type;

template <template <typename> typename Predicate,
          template <typename...> typename Result>
struct filter<Predicate, Result> {
  using type = Result<>;
};

template <template <typename> typename Predicate,
          template <typename...> typename Result, typename T, typename... Ts>
struct filter<Predicate, Result, T, Ts...> {
  template <typename, typename>
  struct append;
  template <typename Head, typename... Tail>
  struct append<Head, Result<Tail...>> {
    using type = Result<Head, Tail...>;
  };

  using type = std::conditional_t<
      Predicate<T>::value,
      typename append<T, typename filter<Predicate, Result, Ts...>::type>::type,
      typename filter<Predicate, Result, Ts...>::type>;
};

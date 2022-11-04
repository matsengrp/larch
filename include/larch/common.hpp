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
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/counted.hpp>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/group_by.hpp>
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
  if constexpr (std::is_same_v<Id, size_t>) {
    if (id >= data.size()) data.resize(id + 1);
    return data[id];
  } else {
    if (id.value >= data.size()) data.resize(id.value + 1);
    return data[id.value];
  }
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
  return ranges::views::transform([&](auto&& i) { return typename DAG::Node{dag, i}; });
}
template <typename DAG>
inline auto ToEdges(DAG dag) {
  return ranges::views::transform([&](auto&& i) { return typename DAG::Edge{dag, i}; });
}

template <typename T>
inline auto To() {
  return ranges::views::transform([&](auto&& i) { return T{i}; });
}
}  // namespace Transform

inline constexpr const auto HashCombine = [](size_t lhs, size_t rhs) noexcept {
  lhs ^= rhs + 0x9e3779b97f4a7c15 + (lhs << 6) + (lhs >> 2);
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

#define MOVE_ONLY(x)           \
  x(x&&) = default;            \
  x(const x&) = delete;        \
  x& operator=(x&&) = default; \
  x& operator=(const x&) = delete

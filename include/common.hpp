#pragma once

#include <cstddef>
#include <limits>
#include <vector>

#include <range/v3/view/transform.hpp>

struct NodeId;
struct EdgeId;
struct CladeIdx;

class DAG;
template <typename>
class NodeView;
using Node = NodeView<const DAG&>;
using MutableNode = NodeView<DAG&>;
template <typename>
class EdgeView;
using Edge = EdgeView<const DAG&>;
using MutableEdge = EdgeView<DAG&>;

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
inline auto ToNodes(const DAG& dag) {
  return ranges::views::transform([&](auto&& i) { return Node{dag, i}; });
}
inline auto ToNodes(DAG& dag) {
  return ranges::views::transform([&](auto&& i) { return MutableNode{dag, i}; });
}
inline auto ToEdges(const DAG& dag) {
  return ranges::views::transform([&](auto&& i) { return Edge{dag, i}; });
}
inline auto ToEdges(DAG& dag) {
  return ranges::views::transform([&](auto&& i) { return MutableEdge{dag, i}; });
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

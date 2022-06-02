#pragma once

#include <cstddef>
#include <limits>
#include <vector>
#include <map>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/subrange.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/all.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/ref.hpp>

static constexpr const size_t NoId = std::numeric_limits<size_t>::max();

struct NodeId {
  size_t value = NoId;
};

struct EdgeId {
  size_t value = NoId;
};

struct CladeIdx {
  size_t value = NoId;
};

struct MutationPosition {
  size_t value = NoId;
};

using Mutations = std::map<MutationPosition, std::pair<char, char>>;

inline bool operator==(NodeId lhs, NodeId rhs) { return lhs.value == rhs.value; }

inline bool operator<(NodeId lhs, NodeId rhs) { return lhs.value < rhs.value; }

inline bool operator==(EdgeId lhs, EdgeId rhs) { return lhs.value == rhs.value; }

inline bool operator<(EdgeId lhs, EdgeId rhs) { return lhs.value < rhs.value; }

inline bool operator==(CladeIdx lhs, CladeIdx rhs) { return lhs.value == rhs.value; }

inline bool operator<(CladeIdx lhs, CladeIdx rhs) { return lhs.value < rhs.value; }

inline bool operator==(MutationPosition lhs, MutationPosition rhs) {
  return lhs.value == rhs.value;
}

inline bool operator<(MutationPosition lhs, MutationPosition rhs) {
  return lhs.value < rhs.value;
}

template <>
struct std::hash<NodeId> {
  size_t operator()(NodeId id) const noexcept { return id.value; }
};

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

class DAG;

template <typename T>
class NodeView;
using Node = NodeView<const DAG&>;
using MutableNode = NodeView<DAG&>;
template <typename T>
class EdgeView;
using Edge = EdgeView<const DAG&>;
using MutableEdge = EdgeView<DAG&>;

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

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

static constexpr size_t NoId = std::numeric_limits<size_t>::max();

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

using Mutations = std::map<MutationPosition, char>;

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

class HistoryDAG;

template <typename T>
class NodeView;
using Node = NodeView<const HistoryDAG&>;
using MutableNode = NodeView<HistoryDAG&>;
template <typename T>
class EdgeView;
using Edge = EdgeView<const HistoryDAG&>;
using MutableEdge = EdgeView<HistoryDAG&>;

namespace Transform {
inline constexpr const auto GetParent = [](auto&& i) { return i.GetParent(); };
inline constexpr const auto GetChild = [](auto&& i) { return i.GetChild(); };
inline constexpr const auto GetId = [](auto&& i) { return i.GetId(); };
}  // namespace Transform

inline constexpr const auto HashCombine = [](size_t lhs, size_t rhs) {
  lhs ^= rhs + 0x9e3779b97f4a7c15 + (lhs << 6) + (lhs >> 2);
  return lhs;
};

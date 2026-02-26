#pragma once

#include <optional>
#include <random>
#include <vector>

#include "larch/dag/dag.hpp"
#include "larch/spr/lca.hpp"

/**
 * @brief Generates random SPR moves on a sampled tree.
 *
 * Given a tree (DAG view), generates random (src, dst, lca) triples that
 * represent valid SPR moves. Uses larch-native NodeIds and DAG traversal
 * — no matOptimize dependency.
 *
 * @tparam DAG The DAG view type (must be a tree)
 */
template <typename DAG>
class RandomMoveGenerator {
 public:
  struct Move {
    NodeId src;
    NodeId dst;
    NodeId lca;
  };

  explicit RandomMoveGenerator(DAG dag,
                               std::optional<uint32_t> seed = std::nullopt);

  /**
   * @brief Generate a random valid SPR move.
   *
   * @param max_distance If >0, only accept moves where the path distance
   *        (path0.size() + path1.size() from LCA) does not exceed this value.
   *        If 0, no distance constraint is applied.
   * @param max_attempts Maximum random attempts before giving up.
   * Returns nullopt if no valid move could be found after max_attempts tries.
   */
  std::optional<Move> GenerateMove(size_t max_distance = 0,
                                   size_t max_attempts = 1000);

 private:
  DAG dag_;
  std::vector<NodeId> searchable_nodes_;  // non-UA nodes
  std::mt19937 gen_;

  bool IsDescendant(NodeId ancestor, NodeId descendant) const;
};

// Implementation

template <typename DAG>
RandomMoveGenerator<DAG>::RandomMoveGenerator(DAG dag,
                                              std::optional<uint32_t> seed)
    : dag_{dag} {
  if (seed.has_value()) {
    gen_.seed(seed.value());
  } else {
    std::random_device rd;
    gen_.seed(rd());
  }

  // Collect all non-UA nodes
  for (auto node : dag_.GetNodes()) {
    if (not node.IsUA()) {
      searchable_nodes_.push_back(node.GetId());
    }
  }
}

template <typename DAG>
bool RandomMoveGenerator<DAG>::IsDescendant(NodeId ancestor,
                                            NodeId descendant) const {
  // Walk from descendant up to root, checking if we hit ancestor
  NodeId current = descendant;
  while (true) {
    if (current == ancestor) {
      return true;
    }
    auto node = dag_.Get(current);
    if (node.IsUA()) {
      return false;
    }
    if (node.GetParentsCount() == 0) {
      return false;
    }
    current = node.GetSingleParent().GetParent().GetId();
  }
}

template <typename DAG>
std::optional<typename RandomMoveGenerator<DAG>::Move>
RandomMoveGenerator<DAG>::GenerateMove(size_t max_distance, size_t max_attempts) {
  if (searchable_nodes_.size() < 3) {
    return std::nullopt;
  }

  std::uniform_int_distribution<size_t> dist(0, searchable_nodes_.size() - 1);

  for (size_t attempt = 0; attempt < max_attempts; attempt++) {
    // Pick random source (must not be tree root)
    NodeId src = searchable_nodes_[dist(gen_)];
    auto src_node = dag_.Get(src);
    if (src_node.IsTreeRoot()) {
      continue;
    }

    // Pick random destination (must not be src, not parent of src,
    // not descendant of src)
    NodeId dst = searchable_nodes_[dist(gen_)];
    if (dst == src) {
      continue;
    }
    auto dst_node = dag_.Get(dst);
    if (dst_node.IsTreeRoot()) {
      continue;
    }

    // dst must not be src's parent
    auto src_parent = src_node.GetSingleParent().GetParent().GetId();
    if (dst == src_parent) {
      continue;
    }

    // dst must not be a descendant of src
    if (IsDescendant(src, dst)) {
      continue;
    }

    // Find LCA
    auto lca_result = FindLCA(dag_.Get(src), dag_.Get(dst));
    NodeId lca = lca_result.lca;

    // Check distance constraint
    if (max_distance > 0) {
      size_t distance = lca_result.path0.size() + lca_result.path1.size();
      if (distance > max_distance) {
        continue;
      }
    }

    // LCA's parent must not be UA (the move must not involve the root edge)
    auto lca_node = dag_.Get(lca);
    if (lca_node.IsUA()) {
      continue;
    }
    if (lca_node.IsTreeRoot()) {
      // LCA is tree root — its parent is UA. This is still a valid move
      // as long as other conditions hold, but we skip for simplicity
      // (matches the mock optimizer behavior).
      continue;
    }

    return Move{src, dst, lca};
  }

  return std::nullopt;
}

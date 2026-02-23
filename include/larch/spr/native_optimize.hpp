#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <functional>
#include <set>

#include "larch/merge/merge.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/spr/spr_view.hpp"
#include "larch/spr/scoring_backend.hpp"
#include "larch/spr/random_optimize.hpp"
#include "larch/parallel/parallel_common.hpp"
#include "larch/parallel/reduction.hpp"
#include "larch/benchmark.hpp"

// ============================================================================
// Data structures for native move enumeration
// ============================================================================

/**
 * @brief DFS traversal indices for O(1) ancestor/descendant checks.
 */
struct DFSInfo {
  size_t dfs_index;
  size_t dfs_end_index;
  size_t level;
};

/**
 * @brief A profitable SPR move found by the native enumerator.
 */
struct ProfitableMove {
  NodeId src;
  NodeId dst;
  NodeId lca;
  int score_change;
};

// ============================================================================
// Fitch helper functions
// ============================================================================

/**
 * @brief Compute the Fitch set from child allele counts using standard
 * Fitch intersection/union semantics.
 */
inline uint8_t FitchSetFromCounts(const std::array<uint8_t, 4>& counts,
                                   uint8_t num_children) {
  if (num_children == 0) return 0;
  uint8_t intersection = 0;
  for (int i = 0; i < 4; i++) {
    if (counts[i] == num_children) {
      intersection |= (1 << i);
    }
  }
  if (intersection) return intersection;
  uint8_t union_set = 0;
  for (int i = 0; i < 4; i++) {
    if (counts[i] > 0) {
      union_set |= (1 << i);
    }
  }
  return union_set;
}

/**
 * @brief Compute Fitch cost at a node from child allele counts.
 * Returns 0 if intersection of children's Fitch sets is non-empty, 1 if empty.
 */
inline int FitchCostFromCounts(const std::array<uint8_t, 4>& counts,
                                uint8_t num_children) {
  if (num_children <= 1) return 0;
  for (int i = 0; i < 4; i++) {
    if (counts[i] == num_children) return 0;
  }
  return 1;
}

// ============================================================================
// TreeIndex: Precomputed tree data for move enumeration
// ============================================================================

/**
 * @brief Precomputed index over a sampled tree for efficient move enumeration.
 *
 * Built once per sampled tree via bottom-up traversal. Stores:
 * - Per-node Fitch sets at variable sites
 * - Per-node child allele counts at variable sites (for incremental updates)
 * - DFS indices for O(1) ancestor checks
 * - List of variable sites and searchable source nodes
 */
template <typename DAG>
class TreeIndex {
 public:
  explicit TreeIndex(DAG dag);

  const std::vector<MutationPosition>& GetVariableSites() const {
    return variable_sites_;
  }
  const std::vector<NodeId>& GetSearchableNodes() const { return searchable_nodes_; }

  bool HasDFSInfo(NodeId node) const {
    return dfs_info_.find(node) != dfs_info_.end();
  }
  const DFSInfo& GetDFSInfo(NodeId node) const { return dfs_info_.at(node); }

  bool IsAncestor(NodeId ancestor, NodeId descendant) const {
    auto anc_it = dfs_info_.find(ancestor);
    auto desc_it = dfs_info_.find(descendant);
    if (anc_it == dfs_info_.end() || desc_it == dfs_info_.end()) return false;
    return anc_it->second.dfs_index <= desc_it->second.dfs_index &&
           desc_it->second.dfs_index < anc_it->second.dfs_end_index;
  }

  /**
   * @brief Get the Fitch set (one-hot bitmask) for a node at a variable site.
   * For leaves, this is the singleton base.
   * For internal nodes, this is the Fitch intersection/union result.
   */
  uint8_t GetFitchSet(NodeId node, size_t site_idx) const {
    auto it = fitch_sets_.find(node);
    if (it == fitch_sets_.end()) return 0;
    return it->second[site_idx];
  }

  /**
   * @brief Get the number of children that have a given allele at a site.
   * counts[allele_bit_index] = number of children with that allele.
   * Only meaningful for internal nodes.
   */
  const std::array<uint8_t, 4>& GetChildCounts(NodeId node, size_t site_idx) const {
    return child_counts_.at(node)[site_idx];
  }

  bool HasChildCounts(NodeId node) const {
    return child_counts_.find(node) != child_counts_.end();
  }

  uint8_t GetNumChildren(NodeId node) const {
    auto it = num_children_.find(node);
    if (it == num_children_.end()) return 0;
    return it->second;
  }

  DAG GetDAG() const { return dag_; }
  NodeId GetTreeRoot() const { return tree_root_; }

  NodeId GetParent(NodeId node) const {
    return dag_.Get(node).GetSingleParent().GetParent().GetId();
  }

  std::vector<NodeId> GetChildren(NodeId node) const {
    std::vector<NodeId> children;
    for (auto clade : dag_.Get(node).GetClades()) {
      for (auto edge : clade) {
        children.push_back(edge.GetChild().GetId());
      }
    }
    return children;
  }

  size_t NumVariableSites() const { return variable_sites_.size(); }

 private:
  void ComputeDFSIndices();
  void DFSVisit(NodeId node, size_t& counter, size_t level);
  void ComputeFitchSets();
  void FitchBottomUp(NodeId node);

  DAG dag_;
  NodeId tree_root_;
  std::vector<MutationPosition> variable_sites_;
  std::vector<NodeId> searchable_nodes_;

  std::unordered_map<NodeId, DFSInfo> dfs_info_;

  // Per-node Fitch sets: fitch_sets_[node][site_idx] = one-hot bitmask
  std::unordered_map<NodeId, std::vector<uint8_t>> fitch_sets_;

  // Per-node child allele counts: child_counts_[node][site_idx][allele] = count
  std::unordered_map<NodeId, std::vector<std::array<uint8_t, 4>>> child_counts_;

  // Per-node actual number of children
  std::unordered_map<NodeId, uint8_t> num_children_;
};

// ============================================================================
// MoveEnumerator: Bounded search for profitable SPR moves
// ============================================================================

/**
 * @brief Exhaustive bounded SPR move enumerator on MADAG.
 *
 * Implements matOptimize's core algorithm: for each source node, walk up
 * to the tree root (bounded by radius), and at each ancestor search the
 * sibling subtrees for improving destinations.
 *
 * The scoring is done by a complete Fitch-based approach: for each candidate
 * move (src, dst), we compute the exact parsimony score change by:
 * 1. Simulating removal of src from its parent
 * 2. Propagating allele count changes upward to the LCA
 * 3. Simulating insertion of src as sibling of dst
 * 4. Propagating those changes upward to the LCA
 * The total score change is the sum of all changes along both paths.
 */
template <typename DAG>
class MoveEnumerator {
 public:
  using Callback = std::function<void(const ProfitableMove&)>;

  explicit MoveEnumerator(const TreeIndex<DAG>& index) : index_{index} {}

  void FindMovesForSource(NodeId src, size_t radius, Callback callback) const;
  void FindAllMoves(size_t radius, Callback callback) const;

  /**
   * @brief Compute the full score change for an SPR move by simulating
   * the Fitch set changes along the removal and insertion paths.
   *
   * Models the actual topology change:
   * - src is pruned from src_parent
   * - If src_parent becomes unifurcation (binary tree), it's collapsed
   * - A new_node with children {src, dst} replaces dst under dst_parent
   * - Fitch sets are recomputed upward on both sides to the LCA
   */
  int ComputeMoveScore(NodeId src, NodeId dst, NodeId lca) const;

 private:

  void UpwardTraversal(NodeId src, size_t radius, Callback& callback) const;

  void SearchSubtree(NodeId node, NodeId src, NodeId lca,
                      size_t radius_left, Callback& callback) const;

  const TreeIndex<DAG>& index_;
};

// ============================================================================
// OptimizeDAGNative: Optimization loop
// ============================================================================

template <typename SampledStorage>
std::vector<RadiusResult> OptimizeDAGNative(Merge& merge,
                                             SampledStorage& sampled_storage,
                                             size_t max_moves_per_radius) {
  auto sampled = sampled_storage.View();
  auto sampled_const = sampled.Const();

  using SampledConstView = decltype(sampled_const);
  using Backend =
      ParsimonyOnlyScoringBackend<std::remove_reference_t<SampledStorage>>;

  std::cout << "  Building tree index..." << std::flush;
  Benchmark index_bench;
  TreeIndex<SampledConstView> tree_index{sampled_const};
  auto index_ms = index_bench.lapMs();
  std::cout << " done (" << index_ms << "ms, "
            << tree_index.GetVariableSites().size() << " variable sites, "
            << tree_index.GetSearchableNodes().size() << " searchable nodes)\n"
            << std::flush;

  size_t max_depth = ComputeTreeMaxDepth(sampled_const);
  size_t max_radius = max_depth * 2;
  std::cout << "  maximum radius is " << max_radius << "\n" << std::flush;

  MoveEnumerator<SampledConstView> enumerator{tree_index};
  std::vector<RadiusResult> results;

  for (size_t rad_exp = 1; (static_cast<size_t>(1) << rad_exp) <= max_radius;
       rad_exp++) {
    size_t radius = static_cast<size_t>(1) << rad_exp;
    std::cout << "  current radius is " << radius << "\n" << std::flush;

    Benchmark radius_bench;

    // Phase 1 (PARALLEL): Enumerate profitable moves
    std::vector<ProfitableMove> all_moves;
    std::mutex moves_mtx;

    auto& searchable = tree_index.GetSearchableNodes();
    ParallelForEach(searchable, [&](NodeId src) {
      enumerator.FindMovesForSource(src, radius,
                                     [&](const ProfitableMove& move) {
                                       std::lock_guard lock{moves_mtx};
                                       all_moves.push_back(move);
                                     });
    });

    // Sort by score_change (most improving first) and cap
    std::sort(all_moves.begin(), all_moves.end(),
              [](const ProfitableMove& a, const ProfitableMove& b) {
                return a.score_change < b.score_change;
              });

    size_t moves_to_apply = std::min(all_moves.size(), max_moves_per_radius);
    all_moves.resize(moves_to_apply);

    auto enumerate_ms = radius_bench.lapMs();

    // Phase 2 (PARALLEL): Create SPR overlays + InitHypotheticalTree
    using SPRStorageType = decltype(AddSPRStorageWithBackend<Backend>(sampled));
    std::vector<std::optional<SPRStorageType>> spr_slots(moves_to_apply);

    std::vector<size_t> work_items(moves_to_apply);
    std::iota(work_items.begin(), work_items.end(), 0);

    ParallelForEach(work_items, [&](size_t idx) {
      auto& move = all_moves[idx];
      spr_slots[idx].emplace(AddSPRStorageWithBackend<Backend>(sampled));
      auto& spr = *spr_slots[idx];
      bool success = spr.View().InitHypotheticalTree(move.src, move.dst, move.lca);
      if (not success) {
        spr_slots[idx].reset();
      }
    });

    auto apply_ms = radius_bench.lapMs();

    // Phase 3 (PARALLEL): Create fragments and overlay MappedNodes
    using FragStorageType =
        decltype(std::declval<SPRStorageType&>().View().MakeFragment());
    std::vector<std::optional<FragStorageType>> frag_slots(moves_to_apply);

    ParallelForEach(work_items, [&](size_t idx) {
      if (not spr_slots[idx].has_value()) {
        return;
      }
      auto& spr = *spr_slots[idx];
      frag_slots[idx].emplace(spr.View().MakeFragment());
      for (auto node : frag_slots[idx]->View().GetNodes()) {
        auto spr_node = spr.View().Get(node.GetId());
        if (not spr_node.IsAppended()) {
          spr_node.template SetOverlay<MappedNodes>();
        }
      }
    });

    auto fragment_ms = radius_bench.lapMs();

    // Phase 4 (BATCH): Merge all fragments
    using FragViewType = decltype(std::declval<FragStorageType&>().View());
    std::vector<FragViewType> frag_views;
    for (size_t idx = 0; idx < moves_to_apply; idx++) {
      if (frag_slots[idx].has_value()) {
        frag_views.push_back(frag_slots[idx]->View());
      }
    }
    if (not frag_views.empty()) {
      merge.AddDAGs(frag_views);
    }
    size_t accepted = frag_views.size();

    auto merge_ms = radius_bench.lapMs();

    std::cout << "  Found " << all_moves.size() << " moves, applied " << accepted
              << ", enumerate=" << enumerate_ms << "ms apply=" << apply_ms
              << "ms frag=" << fragment_ms << "ms merge=" << merge_ms << "ms\n"
              << std::flush;

    results.push_back(RadiusResult{accepted, apply_ms, fragment_ms, merge_ms});
  }

  return results;
}

// ============================================================================
// TreeIndex implementation
// ============================================================================

template <typename DAG>
TreeIndex<DAG>::TreeIndex(DAG dag) : dag_{dag} {
  auto ua = dag_.GetRoot();
  for (auto clade : ua.GetClades()) {
    for (auto edge : clade) {
      tree_root_ = edge.GetChild().GetId();
      break;
    }
    break;
  }

  // Collect variable sites from edge mutations
  std::set<size_t> var_sites_set;
  for (auto edge : dag_.GetEdges()) {
    for (auto& [pos, mut] : edge.GetEdgeMutations()) {
      var_sites_set.insert(pos.value);
    }
  }
  for (auto site : var_sites_set) {
    variable_sites_.push_back({site});
  }

  // Collect searchable nodes (non-UA, non-tree-root)
  for (auto node : dag_.GetNodes()) {
    if (not node.IsUA() and not node.IsTreeRoot()) {
      searchable_nodes_.push_back(node.GetId());
    }
  }

  ComputeDFSIndices();
  ComputeFitchSets();
}

template <typename DAG>
void TreeIndex<DAG>::DFSVisit(NodeId node, size_t& counter, size_t level) {
  DFSInfo info;
  info.dfs_index = counter++;
  info.level = level;

  auto node_view = dag_.Get(node);
  if (not node_view.IsLeaf()) {
    for (auto clade : node_view.GetClades()) {
      for (auto edge : clade) {
        DFSVisit(edge.GetChild().GetId(), counter, level + 1);
      }
    }
  }

  info.dfs_end_index = counter;
  dfs_info_[node] = info;
}

template <typename DAG>
void TreeIndex<DAG>::ComputeDFSIndices() {
  size_t counter = 0;
  DFSVisit(tree_root_, counter, 0);
}

template <typename DAG>
void TreeIndex<DAG>::FitchBottomUp(NodeId node_id) {
  auto node = dag_.Get(node_id);
  size_t n_sites = variable_sites_.size();

  auto& fitch = fitch_sets_[node_id];
  fitch.resize(n_sites);

  if (node.IsLeaf()) {
    for (size_t i = 0; i < n_sites; i++) {
      auto base = node.GetCompactGenome().GetBase(variable_sites_[i],
                                                   dag_.GetReferenceSequence());
      fitch[i] = static_cast<uint8_t>(base_to_singleton(base));
    }
    // No child counts for leaves
  } else {
    std::vector<NodeId> children;
    for (auto clade : node.GetClades()) {
      for (auto edge : clade) {
        auto child_id = edge.GetChild().GetId();
        FitchBottomUp(child_id);
        children.push_back(child_id);
      }
    }

    auto& counts = child_counts_[node_id];
    counts.resize(n_sites);
    num_children_[node_id] = static_cast<uint8_t>(children.size());

    for (size_t i = 0; i < n_sites; i++) {
      auto& c = counts[i];
      c = {0, 0, 0, 0};

      // Count how many children have each allele in their Fitch set
      for (auto child_id : children) {
        uint8_t child_fitch = fitch_sets_[child_id][i];
        for (int j = 0; j < 4; j++) {
          if (child_fitch & (1 << j)) {
            c[j]++;
          }
        }
      }

      fitch[i] = ::FitchSetFromCounts(c, static_cast<uint8_t>(children.size()));
    }
  }
}

template <typename DAG>
void TreeIndex<DAG>::ComputeFitchSets() {
  FitchBottomUp(tree_root_);
}

// ============================================================================
// MoveEnumerator implementation
// ============================================================================

template <typename DAG>
int MoveEnumerator<DAG>::ComputeMoveScore(NodeId src, NodeId dst,
                                           NodeId lca) const {
  size_t n_sites = index_.NumVariableSites();
  if (n_sites == 0) return 0;

  NodeId src_parent = index_.GetParent(src);
  NodeId dst_parent = index_.GetParent(dst);
  bool src_parent_is_binary = (index_.GetNumChildren(src_parent) == 2);

  // Find src's sibling (for unifurcation collapse case)
  NodeId src_sibling{};
  if (src_parent_is_binary) {
    for (auto child : index_.GetChildren(src_parent)) {
      if (child != src) {
        src_sibling = child;
        break;
      }
    }
  }

  // Collect nodes on path from src_parent to LCA (inclusive)
  std::vector<NodeId> src_path;  // src_parent, ..., lca
  {
    NodeId cur = src_parent;
    while (true) {
      src_path.push_back(cur);
      if (cur == lca) break;
      cur = index_.GetParent(cur);
    }
  }

  // Collect nodes on path from dst_parent to LCA (inclusive)
  std::vector<NodeId> dst_path;  // dst_parent, ..., lca
  {
    NodeId cur = dst_parent;
    while (true) {
      dst_path.push_back(cur);
      if (cur == lca) break;
      cur = index_.GetParent(cur);
    }
  }

  // Build a combined list of unique affected nodes, processed bottom-up
  // (deepest first). Use DFS level for ordering.
  std::unordered_set<NodeId> affected_set;
  for (auto n : src_path) affected_set.insert(n);
  for (auto n : dst_path) affected_set.insert(n);

  std::vector<NodeId> affected(affected_set.begin(), affected_set.end());
  std::sort(affected.begin(), affected.end(), [&](NodeId a, NodeId b) {
    return index_.GetDFSInfo(a).level > index_.GetDFSInfo(b).level;
  });

  // For quick path membership checks
  std::unordered_set<NodeId> on_src_path(src_path.begin(), src_path.end());
  std::unordered_set<NodeId> on_dst_path(dst_path.begin(), dst_path.end());

  int total_score_change = 0;

  for (size_t si = 0; si < n_sites; si++) {
    uint8_t src_fitch = index_.GetFitchSet(src, si);
    uint8_t dst_fitch = index_.GetFitchSet(dst, si);

    // New node's Fitch set (children: src, dst)
    uint8_t n_inter = src_fitch & dst_fitch;
    uint8_t new_node_fitch = n_inter ? n_inter : (src_fitch | dst_fitch);
    int new_node_cost = n_inter ? 0 : 1;
    total_score_change += new_node_cost;

    // If src_parent is binary and being collapsed, subtract its old cost
    if (src_parent_is_binary && index_.HasChildCounts(src_parent)) {
      auto counts = index_.GetChildCounts(src_parent, si);
      uint8_t nc = index_.GetNumChildren(src_parent);
      total_score_change -= FitchCostFromCounts(counts, nc);
    }

    // Track new Fitch sets as we process bottom-up
    std::unordered_map<NodeId, uint8_t> new_fitch_map;

    // Seed: if src_parent is binary, it's collapsed to src_sibling's Fitch set
    if (src_parent_is_binary) {
      new_fitch_map[src_parent] = index_.GetFitchSet(src_sibling, si);
    }

    // Process affected nodes bottom-up
    for (auto node : affected) {
      if (src_parent_is_binary && node == src_parent) {
        // Already handled above (collapsed, cost subtracted)
        continue;
      }
      if (not index_.HasChildCounts(node)) continue;

      // Start from original counts
      auto counts = index_.GetChildCounts(node, si);
      uint8_t nc = index_.GetNumChildren(node);
      int old_cost = FitchCostFromCounts(counts, nc);
      uint8_t new_nc = nc;

      // Apply src-side modification to this node's counts
      if (on_src_path.count(node)) {
        if (node == src_parent && !src_parent_is_binary) {
          // Remove src's contribution, decrement child count
          for (int j = 0; j < 4; j++) {
            if (src_fitch & (1 << j)) {
              if (counts[j] > 0) counts[j]--;
            }
          }
          new_nc--;
        } else {
          // Find which child of node is on the src path and has a modified
          // Fitch set. This child is the one closer to src_parent in src_path.
          NodeId src_side_child{};
          for (size_t pi = 0; pi + 1 < src_path.size(); pi++) {
            if (src_path[pi + 1] == node) {
              src_side_child = src_path[pi];
              break;
            }
          }
          if (src_side_child != NodeId{}) {
            auto it = new_fitch_map.find(src_side_child);
            if (it != new_fitch_map.end()) {
              uint8_t old_f = index_.GetFitchSet(src_side_child, si);
              uint8_t new_f = it->second;
              for (int j = 0; j < 4; j++) {
                if (old_f & (1 << j)) {
                  if (counts[j] > 0) counts[j]--;
                }
                if (new_f & (1 << j)) {
                  counts[j]++;
                }
              }
            }
          }
        }
      }

      // Apply dst-side modification to this node's counts
      if (on_dst_path.count(node)) {
        if (node == dst_parent) {
          // Replace dst's contribution with new_node's
          for (int j = 0; j < 4; j++) {
            if (dst_fitch & (1 << j)) {
              if (counts[j] > 0) counts[j]--;
            }
            if (new_node_fitch & (1 << j)) {
              counts[j]++;
            }
          }
        } else {
          // Find which child of node is on the dst path
          NodeId dst_side_child{};
          for (size_t pi = 0; pi + 1 < dst_path.size(); pi++) {
            if (dst_path[pi + 1] == node) {
              dst_side_child = dst_path[pi];
              break;
            }
          }
          if (dst_side_child != NodeId{}) {
            auto it = new_fitch_map.find(dst_side_child);
            if (it != new_fitch_map.end()) {
              uint8_t old_f = index_.GetFitchSet(dst_side_child, si);
              uint8_t new_f = it->second;
              for (int j = 0; j < 4; j++) {
                if (old_f & (1 << j)) {
                  if (counts[j] > 0) counts[j]--;
                }
                if (new_f & (1 << j)) {
                  counts[j]++;
                }
              }
            }
          }
        }
      }

      int new_cost = FitchCostFromCounts(counts, new_nc);
      uint8_t new_fitch = FitchSetFromCounts(counts, new_nc);
      total_score_change += new_cost - old_cost;
      new_fitch_map[node] = new_fitch;
    }
  }

  return total_score_change;
}

template <typename DAG>
void MoveEnumerator<DAG>::SearchSubtree(NodeId node, NodeId src, NodeId lca,
                                         size_t radius_left,
                                         Callback& callback) const {
  if (index_.IsAncestor(src, node) || node == src) {
    return;
  }

  // Compute full move score for (src, node, lca)
  int score = ComputeMoveScore(src, node, lca);
  if (score < 0) {
    callback(ProfitableMove{src, node, lca, score});
  }

  if (radius_left > 0) {
    auto children = index_.GetChildren(node);
    for (auto child : children) {
      SearchSubtree(child, src, lca, radius_left - 1, callback);
    }
  }
}

template <typename DAG>
void MoveEnumerator<DAG>::UpwardTraversal(NodeId src, size_t radius,
                                           Callback& callback) const {
  auto dag = index_.GetDAG();
  auto src_node = dag.Get(src);
  if (src_node.IsTreeRoot() || src_node.IsUA()) {
    return;
  }

  NodeId current = index_.GetParent(src);
  NodeId prev = src;
  size_t levels_up = 0;

  while (levels_up < radius) {
    levels_up++;
    auto current_node = dag.Get(current);

    // Search sibling subtrees
    auto children = index_.GetChildren(current);
    for (auto child : children) {
      if (child == prev) continue;

      size_t remaining_radius = radius - levels_up;
      // For siblings at this level, the LCA is `current`
      SearchSubtree(child, src, current, remaining_radius, callback);
    }

    // Move up
    if (current_node.IsTreeRoot() || current_node.IsUA()) {
      break;
    }
    prev = current;
    current = index_.GetParent(current);
  }
}

template <typename DAG>
void MoveEnumerator<DAG>::FindMovesForSource(NodeId src, size_t radius,
                                              Callback callback) const {
  UpwardTraversal(src, radius, callback);
}

template <typename DAG>
void MoveEnumerator<DAG>::FindAllMoves(size_t radius, Callback callback) const {
  auto& searchable = index_.GetSearchableNodes();
  for (auto src : searchable) {
    FindMovesForSource(src, radius, callback);
  }
}

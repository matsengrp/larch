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

/**
 * @brief Cached result of removing src from its parent and propagating
 * upward to the current LCA level. Reused across all destinations at
 * the same LCA level.
 */
struct SrcRemovalResult {
  int score_change;
  // new_fitch[si]: the new Fitch set of the LCA's child on the src side
  // after src removal has been propagated up
  std::vector<uint8_t> new_fitch;
  // old_fitch[si]: the original Fitch set of that node (for reverting/extending)
  std::vector<uint8_t> old_fitch;
  // Adjustment to the LCA's child count. -1 when src is directly removed
  // from the LCA (level 1), 0 otherwise.
  int lca_nc_adjustment{0};
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
 * - Per-node Fitch sets at variable sites (flat array)
 * - Per-node child allele counts at variable sites (flat array)
 * - Per-node subtree allele unions (flat array, for lower-bound pruning)
 * - DFS indices for O(1) ancestor checks
 * - Precomputed parent/children topology
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
    return node.value < num_nodes_ && has_dfs_info_[node.value];
  }
  const DFSInfo& GetDFSInfo(NodeId node) const { return dfs_info_[node.value]; }

  bool IsAncestor(NodeId ancestor, NodeId descendant) const {
    if (ancestor.value >= num_nodes_ || descendant.value >= num_nodes_) return false;
    if (!has_dfs_info_[ancestor.value] || !has_dfs_info_[descendant.value])
      return false;
    return dfs_info_[ancestor.value].dfs_index <=
               dfs_info_[descendant.value].dfs_index &&
           dfs_info_[descendant.value].dfs_index <
               dfs_info_[ancestor.value].dfs_end_index;
  }

  uint8_t GetFitchSet(NodeId node, size_t site_idx) const {
    return fitch_sets_[node.value * num_variable_sites_ + site_idx];
  }

  const std::array<uint8_t, 4>& GetChildCounts(NodeId node,
                                                 size_t site_idx) const {
    return child_counts_[node.value * num_variable_sites_ + site_idx];
  }

  bool HasChildCounts(NodeId node) const {
    return node.value < num_nodes_ && has_child_counts_[node.value];
  }

  uint8_t GetNumChildren(NodeId node) const { return num_children_[node.value]; }

  uint8_t GetAlleleUnion(NodeId node, size_t site_idx) const {
    return allele_union_[node.value * num_variable_sites_ + site_idx];
  }

  DAG GetDAG() const { return dag_; }
  NodeId GetTreeRoot() const { return tree_root_; }

  NodeId GetParent(NodeId node) const { return parent_[node.value]; }

  const std::vector<NodeId>& GetChildren(NodeId node) const {
    return children_[node.value];
  }

  size_t NumVariableSites() const { return num_variable_sites_; }

 private:
  void ComputeDFSIndices();
  void DFSVisit(NodeId node, size_t& counter, size_t level);
  void ComputeFitchSets();
  void FitchBottomUp(NodeId node);

  DAG dag_;
  NodeId tree_root_;
  std::vector<MutationPosition> variable_sites_;
  std::vector<NodeId> searchable_nodes_;

  size_t num_nodes_{0};
  size_t num_variable_sites_{0};

  // Flat arrays indexed by node.value
  std::vector<DFSInfo> dfs_info_;
  std::vector<bool> has_dfs_info_;

  // fitch_sets_[node.value * num_variable_sites_ + site_idx]
  std::vector<uint8_t> fitch_sets_;

  // child_counts_[node.value * num_variable_sites_ + site_idx][allele]
  std::vector<std::array<uint8_t, 4>> child_counts_;
  std::vector<bool> has_child_counts_;

  std::vector<uint8_t> num_children_;

  // allele_union_[node.value * num_variable_sites_ + site_idx]:
  // OR of all alleles in this node's subtree at this site
  std::vector<uint8_t> allele_union_;

  // Precomputed topology
  std::vector<NodeId> parent_;
  std::vector<std::vector<NodeId>> children_;
};

// ============================================================================
// MoveEnumerator: Bounded search for profitable SPR moves
// ============================================================================

/**
 * @brief Exhaustive bounded SPR move enumerator on MADAG with incremental
 * scoring and lower-bound pruning.
 *
 * For each source node, walks up toward the tree root (bounded by radius).
 * At each ancestor level, caches the src-side removal state and searches
 * sibling subtrees, pruning branches where a lower bound on the score
 * change exceeds the best score found so far.
 */
template <typename DAG>
class MoveEnumerator {
 public:
  using Callback = std::function<void(const ProfitableMove&)>;

  explicit MoveEnumerator(const TreeIndex<DAG>& index) : index_{index} {}

  void FindMovesForSource(NodeId src, size_t radius, Callback callback) const;
  void FindAllMoves(size_t radius, Callback callback) const;

  /**
   * @brief Compute score change independently (no caching). Kept for
   * verification in tests.
   */
  int ComputeMoveScore(NodeId src, NodeId dst, NodeId lca) const;

  SrcRemovalResult ComputeInitialRemoval(NodeId src) const;

  void PropagateRemovalUpward(NodeId current_lca,
                               SrcRemovalResult& removal) const;

  int ComputeLowerBound(NodeId src, const SrcRemovalResult& removal,
                         NodeId subtree_root) const;

  int ComputeMoveScoreCached(NodeId src, NodeId dst, NodeId lca,
                              const SrcRemovalResult& removal) const;

 private:
  void UpwardTraversal(NodeId src, size_t radius, Callback& callback) const;

  void SearchSubtreeWithBound(NodeId node, NodeId src, NodeId lca,
                               size_t radius_left,
                               const SrcRemovalResult& removal,
                               int& best_score, Callback& callback) const;

  void SearchSubtreeDirectly(NodeId node, NodeId src, NodeId lca,
                              size_t radius_left, int& best_score,
                              Callback& callback) const;

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
  num_variable_sites_ = variable_sites_.size();

  // Determine num_nodes_ from DAG node count
  num_nodes_ = dag_.GetNodesCount();

  // Allocate flat arrays
  dfs_info_.resize(num_nodes_);
  has_dfs_info_.assign(num_nodes_, false);

  fitch_sets_.resize(num_nodes_ * num_variable_sites_, 0);
  child_counts_.resize(num_nodes_ * num_variable_sites_, {0, 0, 0, 0});
  has_child_counts_.assign(num_nodes_, false);
  num_children_.assign(num_nodes_, 0);
  allele_union_.resize(num_nodes_ * num_variable_sites_, 0);

  parent_.resize(num_nodes_);
  children_.resize(num_nodes_);

  // Build precomputed topology
  for (auto node : dag_.GetNodes()) {
    NodeId nid = node.GetId();
    if (node.IsUA()) continue;

    // Parent
    if (!node.IsTreeRoot()) {
      parent_[nid.value] = node.GetSingleParent().GetParent().GetId();
    }

    // Children
    for (auto clade : node.GetClades()) {
      for (auto edge : clade) {
        children_[nid.value].push_back(edge.GetChild().GetId());
      }
    }
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

  for (auto child : children_[node.value]) {
    DFSVisit(child, counter, level + 1);
  }

  info.dfs_end_index = counter;
  dfs_info_[node.value] = info;
  has_dfs_info_[node.value] = true;
}

template <typename DAG>
void TreeIndex<DAG>::ComputeDFSIndices() {
  size_t counter = 0;
  DFSVisit(tree_root_, counter, 0);
}

template <typename DAG>
void TreeIndex<DAG>::FitchBottomUp(NodeId node_id) {
  auto node = dag_.Get(node_id);
  size_t base_offset = node_id.value * num_variable_sites_;

  if (node.IsLeaf()) {
    for (size_t i = 0; i < num_variable_sites_; i++) {
      auto base = node.GetCompactGenome().GetBase(variable_sites_[i],
                                                   dag_.GetReferenceSequence());
      uint8_t singleton = static_cast<uint8_t>(base_to_singleton(base));
      fitch_sets_[base_offset + i] = singleton;
      allele_union_[base_offset + i] = singleton;
    }
  } else {
    auto& children = children_[node_id.value];
    for (auto child_id : children) {
      FitchBottomUp(child_id);
    }

    has_child_counts_[node_id.value] = true;
    uint8_t nc = static_cast<uint8_t>(children.size());
    num_children_[node_id.value] = nc;

    for (size_t i = 0; i < num_variable_sites_; i++) {
      auto& c = child_counts_[base_offset + i];
      c = {0, 0, 0, 0};

      uint8_t au = 0;
      for (auto child_id : children) {
        size_t child_offset = child_id.value * num_variable_sites_ + i;
        uint8_t child_fitch = fitch_sets_[child_offset];
        for (int j = 0; j < 4; j++) {
          if (child_fitch & (1 << j)) {
            c[j]++;
          }
        }
        au |= allele_union_[child_offset];
      }

      fitch_sets_[base_offset + i] = ::FitchSetFromCounts(c, nc);
      allele_union_[base_offset + i] = au;
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

  NodeId src_sibling{};
  if (src_parent_is_binary) {
    for (auto child : index_.GetChildren(src_parent)) {
      if (child != src) {
        src_sibling = child;
        break;
      }
    }
  }

  std::vector<NodeId> src_path;
  {
    NodeId cur = src_parent;
    while (true) {
      src_path.push_back(cur);
      if (cur == lca) break;
      cur = index_.GetParent(cur);
    }
  }

  std::vector<NodeId> dst_path;
  {
    NodeId cur = dst_parent;
    while (true) {
      dst_path.push_back(cur);
      if (cur == lca) break;
      cur = index_.GetParent(cur);
    }
  }

  std::unordered_set<NodeId> affected_set;
  for (auto n : src_path) affected_set.insert(n);
  for (auto n : dst_path) affected_set.insert(n);

  std::vector<NodeId> affected(affected_set.begin(), affected_set.end());
  std::sort(affected.begin(), affected.end(), [&](NodeId a, NodeId b) {
    return index_.GetDFSInfo(a).level > index_.GetDFSInfo(b).level;
  });

  std::unordered_set<NodeId> on_src_path(src_path.begin(), src_path.end());
  std::unordered_set<NodeId> on_dst_path(dst_path.begin(), dst_path.end());

  int total_score_change = 0;

  for (size_t si = 0; si < n_sites; si++) {
    uint8_t src_fitch = index_.GetFitchSet(src, si);
    uint8_t dst_fitch = index_.GetFitchSet(dst, si);

    uint8_t n_inter = src_fitch & dst_fitch;
    uint8_t new_node_fitch = n_inter ? n_inter : (src_fitch | dst_fitch);
    int new_node_cost = n_inter ? 0 : 1;
    total_score_change += new_node_cost;

    if (src_parent_is_binary && index_.HasChildCounts(src_parent)) {
      auto counts = index_.GetChildCounts(src_parent, si);
      uint8_t nc = index_.GetNumChildren(src_parent);
      total_score_change -= FitchCostFromCounts(counts, nc);
    }

    std::unordered_map<NodeId, uint8_t> new_fitch_map;

    if (src_parent_is_binary) {
      new_fitch_map[src_parent] = index_.GetFitchSet(src_sibling, si);
    }

    for (auto node : affected) {
      if (src_parent_is_binary && node == src_parent) {
        continue;
      }
      if (not index_.HasChildCounts(node)) continue;

      auto counts = index_.GetChildCounts(node, si);
      uint8_t nc = index_.GetNumChildren(node);
      int old_cost = FitchCostFromCounts(counts, nc);
      uint8_t new_nc = nc;

      if (on_src_path.count(node)) {
        if (node == src_parent && !src_parent_is_binary) {
          for (int j = 0; j < 4; j++) {
            if (src_fitch & (1 << j)) {
              if (counts[j] > 0) counts[j]--;
            }
          }
          new_nc--;
        } else {
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

      if (on_dst_path.count(node)) {
        if (node == dst_parent) {
          for (int j = 0; j < 4; j++) {
            if (dst_fitch & (1 << j)) {
              if (counts[j] > 0) counts[j]--;
            }
            if (new_node_fitch & (1 << j)) {
              counts[j]++;
            }
          }
        } else {
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

// ============================================================================
// Incremental removal + lower-bound pruning
// ============================================================================

template <typename DAG>
SrcRemovalResult MoveEnumerator<DAG>::ComputeInitialRemoval(NodeId src) const {
  size_t n_sites = index_.NumVariableSites();
  SrcRemovalResult result;
  result.score_change = 0;
  result.new_fitch.resize(n_sites);
  result.old_fitch.resize(n_sites);

  NodeId src_parent = index_.GetParent(src);
  bool src_parent_is_binary = (index_.GetNumChildren(src_parent) == 2);

  if (src_parent_is_binary) {
    // src_parent gets collapsed: find sibling
    NodeId sibling{};
    for (auto child : index_.GetChildren(src_parent)) {
      if (child != src) {
        sibling = child;
        break;
      }
    }

    for (size_t si = 0; si < n_sites; si++) {
      // Subtract old src_parent cost
      auto& counts = index_.GetChildCounts(src_parent, si);
      uint8_t nc = index_.GetNumChildren(src_parent);
      result.score_change -= FitchCostFromCounts(counts, nc);

      // After collapse, src_parent's Fitch set becomes sibling's
      result.old_fitch[si] = index_.GetFitchSet(src_parent, si);
      result.new_fitch[si] = index_.GetFitchSet(sibling, si);
    }
  } else {
    // src_parent loses one child
    for (size_t si = 0; si < n_sites; si++) {
      auto counts = index_.GetChildCounts(src_parent, si);
      uint8_t nc = index_.GetNumChildren(src_parent);
      int old_cost = FitchCostFromCounts(counts, nc);

      // Remove src's contribution
      uint8_t src_fitch = index_.GetFitchSet(src, si);
      for (int j = 0; j < 4; j++) {
        if (src_fitch & (1 << j)) {
          if (counts[j] > 0) counts[j]--;
        }
      }
      uint8_t new_nc = nc - 1;
      int new_cost = FitchCostFromCounts(counts, new_nc);
      uint8_t new_fitch = FitchSetFromCounts(counts, new_nc);

      result.score_change += new_cost - old_cost;
      result.old_fitch[si] = index_.GetFitchSet(src_parent, si);
      result.new_fitch[si] = new_fitch;
    }
  }

  return result;
}

template <typename DAG>
void MoveEnumerator<DAG>::PropagateRemovalUpward(
    NodeId current_lca, SrcRemovalResult& removal) const {
  // current_lca is the node we're propagating through.
  // removal.new_fitch/old_fitch currently describe the child of current_lca
  // on the src side. We need to update current_lca's counts and then
  // set new_fitch/old_fitch to describe current_lca itself.
  size_t n_sites = index_.NumVariableSites();

  for (size_t si = 0; si < n_sites; si++) {
    auto counts = index_.GetChildCounts(current_lca, si);
    uint8_t nc = index_.GetNumChildren(current_lca);
    int old_cost = FitchCostFromCounts(counts, nc);

    uint8_t old_child_f = removal.old_fitch[si];
    uint8_t new_child_f = removal.new_fitch[si];

    // Update counts: remove old child's contribution, add new
    for (int j = 0; j < 4; j++) {
      if (old_child_f & (1 << j)) {
        if (counts[j] > 0) counts[j]--;
      }
      if (new_child_f & (1 << j)) {
        counts[j]++;
      }
    }

    int new_cost = FitchCostFromCounts(counts, nc);
    uint8_t new_fitch = FitchSetFromCounts(counts, nc);

    removal.score_change += new_cost - old_cost;
    removal.old_fitch[si] = index_.GetFitchSet(current_lca, si);
    removal.new_fitch[si] = new_fitch;
  }
}

template <typename DAG>
int MoveEnumerator<DAG>::ComputeLowerBound(
    NodeId src, const SrcRemovalResult& /*removal*/,
    NodeId subtree_root) const {
  size_t n_sites = index_.NumVariableSites();
  // Count sites where src's allele doesn't appear in the subtree.
  // At these sites, for ANY destination dst in the subtree,
  // src_fitch & dst_fitch == 0, so new_node costs +1.
  // This is a valid lower bound on the new_node cost component.
  // We don't include removal.score_change because dst-side and LCA
  // effects can compensate for it, making total score < lower_bound.
  int forced_new_node_cost = 0;
  for (size_t si = 0; si < n_sites; si++) {
    uint8_t src_fitch = index_.GetFitchSet(src, si);
    uint8_t subtree_alleles = index_.GetAlleleUnion(subtree_root, si);
    if (!(src_fitch & subtree_alleles)) {
      forced_new_node_cost += 1;
    }
  }
  return forced_new_node_cost;
}

template <typename DAG>
int MoveEnumerator<DAG>::ComputeMoveScoreCached(
    NodeId src, NodeId dst, NodeId lca,
    const SrcRemovalResult& removal) const {
  size_t n_sites = index_.NumVariableSites();
  if (n_sites == 0) return 0;

  NodeId dst_parent = index_.GetParent(dst);

  // Collect dst path from dst_parent to LCA (inclusive)
  std::vector<NodeId> dst_path;
  {
    NodeId cur = dst_parent;
    while (true) {
      dst_path.push_back(cur);
      if (cur == lca) break;
      cur = index_.GetParent(cur);
    }
  }

  int total = removal.score_change;

  for (size_t si = 0; si < n_sites; si++) {
    uint8_t src_fitch = index_.GetFitchSet(src, si);
    uint8_t dst_fitch = index_.GetFitchSet(dst, si);

    // New node cost (children: src, dst)
    uint8_t n_inter = src_fitch & dst_fitch;
    uint8_t new_node_fitch = n_inter ? n_inter : (src_fitch | dst_fitch);
    int new_node_cost = n_inter ? 0 : 1;
    total += new_node_cost;

    // Process dst path bottom-up, tracking new Fitch sets
    uint8_t prev_new_fitch = new_node_fitch;
    uint8_t prev_old_fitch = dst_fitch;

    for (size_t pi = 0; pi < dst_path.size(); pi++) {
      NodeId node = dst_path[pi];
      if (!index_.HasChildCounts(node)) continue;

      auto counts = index_.GetChildCounts(node, si);
      uint8_t nc = index_.GetNumChildren(node);
      int old_cost = FitchCostFromCounts(counts, nc);

      // At dst_parent: replace dst with new_node
      // At higher nodes: replace the modified child's old fitch with new fitch
      uint8_t old_child_f;
      uint8_t new_child_f;
      if (pi == 0) {
        old_child_f = dst_fitch;
        new_child_f = new_node_fitch;
      } else {
        old_child_f = prev_old_fitch;
        new_child_f = prev_new_fitch;
      }

      for (int j = 0; j < 4; j++) {
        if (old_child_f & (1 << j)) {
          if (counts[j] > 0) counts[j]--;
        }
        if (new_child_f & (1 << j)) {
          counts[j]++;
        }
      }

      // At LCA: also apply the src-side removal effect
      uint8_t effective_nc = nc;
      if (node == lca) {
        uint8_t src_old_f = removal.old_fitch[si];
        uint8_t src_new_f = removal.new_fitch[si];
        for (int j = 0; j < 4; j++) {
          if (src_old_f & (1 << j)) {
            if (counts[j] > 0) counts[j]--;
          }
          if (src_new_f & (1 << j)) {
            counts[j]++;
          }
        }
        effective_nc = static_cast<uint8_t>(
            static_cast<int>(nc) + removal.lca_nc_adjustment);
      }

      int new_cost = FitchCostFromCounts(counts, effective_nc);
      uint8_t new_fitch = FitchSetFromCounts(counts, effective_nc);
      total += new_cost - old_cost;

      prev_old_fitch = index_.GetFitchSet(node, si);
      prev_new_fitch = new_fitch;
    }
  }

  // Subtract the src-side path costs that were already counted in
  // removal.score_change for nodes between src_parent and LCA.
  // But we need to undo the LCA portion of removal.score_change since
  // we've recomputed it above with both src and dst effects combined.
  //
  // Actually, removal.score_change includes cost changes from src_parent
  // up to (but NOT including) the LCA. The LCA is the topmost node in
  // the traversal, and PropagateRemovalUpward was called for nodes
  // *below* LCA. So we need to check what's in removal.
  //
  // Wait — let's reconsider. The removal state when UpwardTraversal calls
  // us has been propagated through all nodes from src_parent up to (but
  // not including) the current LCA. The new_fitch/old_fitch describe
  // the LCA's child on the src side. So removal.score_change does NOT
  // include LCA's cost change. The dst path loop above processes LCA and
  // applies both src-side and dst-side effects there. This is correct.

  return total;
}

template <typename DAG>
void MoveEnumerator<DAG>::SearchSubtreeWithBound(
    NodeId node, NodeId src, NodeId lca, size_t radius_left,
    const SrcRemovalResult& removal, int& best_score,
    Callback& callback) const {
  if (index_.IsAncestor(src, node) || node == src) {
    return;
  }

  // Score this destination
  int score = ComputeMoveScoreCached(src, node, lca, removal);
  if (score < 0) {
    callback(ProfitableMove{src, node, lca, score});
    if (score < best_score) {
      best_score = score;
    }
  }

  if (radius_left == 0) return;

  for (auto child : index_.GetChildren(node)) {
    SearchSubtreeWithBound(child, src, lca, radius_left - 1, removal, best_score,
                            callback);
  }
}

template <typename DAG>
void MoveEnumerator<DAG>::SearchSubtreeDirectly(
    NodeId node, NodeId src, NodeId lca, size_t radius_left,
    int& best_score, Callback& callback) const {
  if (index_.IsAncestor(src, node) || node == src) {
    return;
  }

  int score = ComputeMoveScore(src, node, lca);
  if (score < 0) {
    callback(ProfitableMove{src, node, lca, score});
    if (score < best_score) {
      best_score = score;
    }
  }

  if (radius_left > 0) {
    for (auto child : index_.GetChildren(node)) {
      SearchSubtreeDirectly(child, src, lca, radius_left - 1, best_score,
                             callback);
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

  size_t n_sites = index_.NumVariableSites();

  // Level-1 removal: src is directly removed from LCA (= src_parent).
  // score_change = 0 (no intermediate nodes), old/new_fitch describe
  // src's contribution, nc_adjustment = -1.
  SrcRemovalResult removal;
  removal.score_change = 0;
  removal.old_fitch.resize(n_sites);
  removal.new_fitch.resize(n_sites, 0);
  removal.lca_nc_adjustment = -1;
  for (size_t si = 0; si < n_sites; si++) {
    removal.old_fitch[si] = index_.GetFitchSet(src, si);
  }

  NodeId current = index_.GetParent(src);
  NodeId prev = src;
  size_t levels_up = 0;
  int best_score = 0;

  while (levels_up < radius) {
    levels_up++;
    auto current_node = dag.Get(current);

    auto& children = index_.GetChildren(current);
    size_t remaining_radius = radius - levels_up;

    for (auto child : children) {
      if (child == prev) continue;
      SearchSubtreeWithBound(child, src, current, remaining_radius, removal,
                              best_score, callback);
    }

    // Move up
    if (current_node.IsTreeRoot() || current_node.IsUA()) {
      break;
    }

    if (levels_up == 1) {
      // Transition from level 1 to level 2.
      // ComputeInitialRemoval gives the state describing src_parent's
      // fitch change as a child of its parent (for level 2+).
      removal = ComputeInitialRemoval(src);
      removal.lca_nc_adjustment = 0;
    } else {
      PropagateRemovalUpward(current, removal);
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

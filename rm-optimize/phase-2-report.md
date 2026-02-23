# Phase 2: Fast Native Move Enumeration — Report

## Overview

Phase 2 optimized the native SPR move enumeration implemented in Phase 1. The primary
changes are flat array indexing (replacing hash maps), precomputed topology, and
incremental cached scoring that avoids redundant computation across destinations sharing
the same source node and LCA.

All changes are in two files:
- `include/larch/spr/native_optimize.hpp`
- `test/test_native_optimize.cpp`

## Benchmark Results (seedtree, debug build, clean system)

Dataset: seedtree (~1040 nodes, 1238 variable sites, 1036 searchable nodes)

### Enumeration time per radius

| Radius | Phase 1 | Phase 2 | Speedup |
|--------|---------|---------|---------|
| 2 | — | 72ms | — |
| 4 | — | 309ms | — |
| 8 | 2,300ms | 1,815ms | 1.3x |
| 16 | 30,400ms | 10,054ms | 3.0x |
| 32 | 128,700ms | 47,873ms | 2.7x |
| 64 | killed | 72,558ms | completes |

### Per-iteration total optimize time

| | Phase 1 | Phase 2 |
|---|---------|---------|
| All radii (2-64) | >300s (killed at r64) | ~163s |

### 3-iteration total wall time

Phase 2: **8m12s** (492s). Highly consistent across iterations (162s, 164s, 164s).

## Implementation Details

### 1. Flat Arrays in TreeIndex

Replaced all `std::unordered_map<NodeId, ...>` members with flat `std::vector` arrays
indexed by `node.value`. Node IDs in sampled trees are contiguous 0..N-1, making this
safe and cache-friendly.

**Before:**
```cpp
std::unordered_map<NodeId, DFSInfo> dfs_info_;
std::unordered_map<NodeId, std::vector<uint8_t>> fitch_sets_;
std::unordered_map<NodeId, std::vector<std::array<uint8_t, 4>>> child_counts_;
std::unordered_map<NodeId, uint8_t> num_children_;
```

**After:**
```cpp
size_t num_nodes_;
size_t num_variable_sites_;
std::vector<DFSInfo> dfs_info_;           // [num_nodes_]
std::vector<bool> has_dfs_info_;          // [num_nodes_]
std::vector<uint8_t> fitch_sets_;         // [num_nodes_ * num_variable_sites_]
std::vector<std::array<uint8_t, 4>> child_counts_;  // [num_nodes_ * num_variable_sites_]
std::vector<bool> has_child_counts_;      // [num_nodes_]
std::vector<uint8_t> num_children_;       // [num_nodes_]
std::vector<uint8_t> allele_union_;       // [num_nodes_ * num_variable_sites_]
std::vector<NodeId> parent_;              // [num_nodes_]
std::vector<std::vector<NodeId>> children_;  // [num_nodes_]
```

Accessors become O(1) flat index lookups:
```cpp
uint8_t GetFitchSet(NodeId node, size_t site_idx) const {
  return fitch_sets_[node.value * num_variable_sites_ + site_idx];
}
```

Memory: ~1040 nodes x ~1238 sites x ~6 bytes/entry = ~8 MB total. Affordable.

`GetParent()` and `GetChildren()` were also converted from DAG traversal to precomputed
flat arrays, eliminating per-call overhead in the hot loop.

### 2. Allele Union (Subtree Allele Sets)

Added `allele_union_[node][site]`: the bitwise OR of all leaf alleles in a node's subtree
at each variable site. Computed during `FitchBottomUp` alongside Fitch sets.

For leaves: `allele_union = base_to_singleton(leaf_base)`
For internal nodes: `allele_union = OR of children's allele_union`

This provides the information needed for lower-bound pruning (see section 5).

### 3. SrcRemovalResult and Incremental Scoring

The key optimization. In Phase 1, every `ComputeMoveScore(src, dst, lca)` independently
rebuilt the full src-to-LCA path. For a source with 500 potential destinations at the
same LCA level, this means 500 redundant path computations.

Phase 2 introduces `SrcRemovalResult`:
```cpp
struct SrcRemovalResult {
  int score_change;              // accumulated cost delta from removing src
  std::vector<uint8_t> new_fitch; // new Fitch set of LCA's src-side child
  std::vector<uint8_t> old_fitch; // original Fitch set (for propagation)
  int lca_nc_adjustment;         // -1 at level 1 (src directly removed), 0 otherwise
};
```

Built once per source via `ComputeInitialRemoval(src)`, then propagated incrementally
via `PropagateRemovalUpward(node, removal)` as the search walks up the tree. The state
is reused for all destinations at each LCA level.

`ComputeMoveScoreCached(src, dst, lca, removal)` then only computes the dst-side path
fresh, taking the src-side from the cached removal state. At the LCA, both effects are
combined.

### 4. Level-1 Cached Scoring

A non-trivial challenge: at level 1, the LCA equals src_parent, meaning src is directly
removed as a child (not replaced by a modified subtree). This changes the child count
at the LCA.

**Discovery:** The `lca_nc_adjustment` field in `SrcRemovalResult` handles this. At
level 1, it's set to -1 (src removed), and `ComputeMoveScoreCached` uses
`effective_nc = nc + lca_nc_adjustment` when computing costs at the LCA node.

**Key bug found during development:** For binary src_parent, an early attempt put
`-FitchCost(src_parent)` into `removal.score_change` AND had `ComputeMoveScoreCached`
compute `new_cost - old_cost` at the LCA, double-counting the old cost subtraction. The
fix: level-1 removal has `score_change = 0` (no intermediate nodes), and the binary
collapse cost is subtracted separately (not in the per-site LCA loop but as a one-time
deduction in `score_change`).

**Transition from level 1 to level 2:** After searching at level 1, the removal state is
replaced by `ComputeInitialRemoval(src)` which describes src_parent's fitch change as
a child of its parent (the level-2 LCA). `lca_nc_adjustment` returns to 0.

### 5. Lower-Bound Pruning (Infrastructure Only)

`ComputeLowerBound(src, removal, subtree_root)` counts sites where src's allele doesn't
appear anywhere in the subtree's `allele_union`. At these sites, the new_node (children:
src, dst) must pay +1 regardless of which dst is chosen.

**Discovery: valid lower bounds are hard.** The forced new_node cost is a correct lower
bound on the new_node cost component, but NOT on the total score change. The total
includes removal path costs, dst-side path costs, and LCA combined effects, all of which
can be negative and compensate for the new_node costs. Attempts to use
`removal.score_change + forced_new_node_cost` as a bound caused incorrect pruning:

- Example: `removal.score_change = 2, forced_new_node_cost = 0, bound = 2`. But actual
  best score in subtree was -1 (profitable move exists). The bound incorrectly said no
  profitable move possible, because the dst-side path saved more than the src-side cost.

- The pruning condition `if (lb >= best_score && best_score < 0) continue` would prune
  subtrees containing moves worse than the current best but still profitable — missing
  moves that the optimization loop could use.

**Current state:** The pruning infrastructure (`allele_union_`, `ComputeLowerBound`) is
in place but pruning is disabled in the search. Enabling it requires either:
(a) A tighter bound that accounts for dst-side compensation, or
(b) Accepting that pruning may miss some profitable moves (only finding the best per
    source, not all profitable ones)

### 6. SearchSubtreeDirectly (Kept for ComputeMoveScore verification)

`SearchSubtreeDirectly` uses the independent `ComputeMoveScore` (no caching). It's kept
as a code path for test verification but no longer used in the main enumeration.

## Test Changes

Two new tests added:

**Test 7: Pruned vs exhaustive** — Enumerates all moves via the pruned/cached path
(`FindAllMoves`) and independently via exhaustive `ComputeMoveScore` on the same tree
structure. Verifies the same set of (src, dst) pairs are found. Currently: 10 = 10,
0 missing, 0 extra.

**Test 8: Cached vs independent scoring** — For every (src, dst, lca) reachable at levels
1 and 2, verifies `ComputeMoveScoreCached` returns the same score as `ComputeMoveScore`.
Currently: 48 checks, 0 mismatches.

## Architecture Decisions

### Why not prune?

The plan called for aggressive pruning using `removal.score_change + forced_new_node_cost`
as a lower bound. This turned out to be incorrect — see section 5. A correct lower bound
on the total score would need to account for the maximum possible savings from dst-side
path changes and LCA combined effects, which is difficult to compute cheaply.

The conservative choice (no pruning) ensures we find all profitable moves. The cached
scoring alone provides 2.7-3x speedup at large radii, and radius 64 now completes.

### Why level-1 cached scoring matters

Level 1 (LCA = src_parent) accounts for the majority of destination nodes at large radii.
A source's siblings span the entire opposite subtree of the tree. Without cached scoring
at level 1, most destinations still use the expensive `ComputeMoveScore`. Adding level-1
support required the `lca_nc_adjustment` mechanism but provides the largest single speedup.

### Memory overhead

The flat arrays add ~8 MB for the seedtree dataset (1040 nodes, 1238 sites). This is
negligible compared to the DAG's own memory usage (~200 MB peak).

## Remaining Bottlenecks

1. **Enumeration still dominates** — At radius 64, enumeration is 72s vs 5s for apply +
   merge. The O(n_src * n_dst_per_radius * n_sites) complexity hasn't changed
   asymptotically, only the constant factor.

2. **No pruning** — The allele_union infrastructure is ready but needs a tighter bound
   to be useful. One approach: accept missing some non-best moves and only guarantee
   finding the best move per source.

3. **Per-destination allocation** — `ComputeMoveScoreCached` allocates a `dst_path`
   vector per call. This could be preallocated per thread.

4. **Binary src_parent handling** — The binary collapse case adds complexity to both
   `ComputeInitialRemoval` and the LCA scoring logic. A refactor to handle unifurcation
   more uniformly could simplify the code.

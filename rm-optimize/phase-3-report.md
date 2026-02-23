# Phase 3: Eliminate Allocations + Micro-Optimizations — Report

## Overview

Phase 3 eliminated per-call heap allocations from the hot path and added micro-optimizations
to `ComputeMoveScoreCached`, the function called once per (src, dst) pair — millions of
times per radius.

All changes are in one file:
- `include/larch/spr/native_optimize.hpp`

No test changes were needed; all existing tests pass without modification.

## Benchmark Results (seedtree, debug build, clean system)

Dataset: seedtree (~1040 nodes, 1238 variable sites, 1036 searchable nodes)

### Enumeration time per radius

| Radius | Phase 2 | Phase 3 | Speedup |
|--------|---------|---------|---------|
| 2 | 72ms | 55ms | 1.31x |
| 4 | 309ms | 228ms | 1.36x |
| 8 | 1,815ms | 1,353ms | 1.34x |
| 16 | 10,054ms | 7,681ms | 1.31x |
| 32 | 47,873ms | 36,811ms | 1.30x |
| 64 | 72,558ms | 56,310ms | 1.29x |

### Per-iteration total enumerate time

| | Phase 2 | Phase 3 |
|---|---------|---------|
| All radii (2-64) | ~132s | ~102s |

### 3-iteration total wall time

| | Phase 2 | Phase 3 |
|---|---------|---------|
| Total | 492s (8m12s) | 401s (6m41s) |

### Cumulative speedup from Phase 1

| Radius | Phase 1 | Phase 3 | Cumulative speedup |
|--------|---------|---------|-------------------|
| 8 | 2,300ms | 1,353ms | 1.7x |
| 16 | 30,400ms | 7,681ms | 4.0x |
| 32 | 128,700ms | 36,811ms | 3.5x |
| 64 | killed | 56,310ms | ∞ (now completes) |

## Implementation Details

### 1. ScratchBuffers (eliminate per-call allocation)

Added a `ScratchBuffers` struct allocated once per `FindMovesForSource` call:

```cpp
struct ScratchBuffers {
  std::vector<uint8_t> new_node_fitch;   // [n_sites]
  std::vector<uint8_t> prev_old_fitch;   // [n_sites]
  std::vector<uint8_t> prev_new_fitch;   // [n_sites]
  SrcRemovalResult removal;              // pre-allocated removal state
};
```

Threaded through the entire call chain: `FindMovesForSource` → `UpwardTraversal` →
`SearchSubtreeWithBound` → `ComputeMoveScoreCached`. Each source node's search reuses
the same buffers for all destination evaluations.

**Impact:** Eliminates millions of small heap allocations per radius. Each
`ComputeMoveScoreCached` call previously allocated a `std::vector<NodeId>` for the dst
path (~10-20 entries). Each level transition in `UpwardTraversal` allocated two 1238-element
`std::vector<uint8_t>` for `SrcRemovalResult`.

### 2. Parent-pointer traversal (eliminate dst_path vector)

Replaced the two-pass approach (build `std::vector<NodeId> dst_path`, then iterate) with
a single-pass parent-pointer walk:

**Before:**
```cpp
std::vector<NodeId> dst_path;
NodeId cur = dst_parent;
while (true) { dst_path.push_back(cur); if (cur == lca) break; cur = GetParent(cur); }
// ...
for (size_t pi = 0; pi < dst_path.size(); pi++) { /* process dst_path[pi] */ }
```

**After:**
```cpp
NodeId node = dst_parent;
bool is_first = true;
while (true) {
  // process node, using prev_old/prev_new from scratch buffers
  if (node == lca) break;
  node = index_.GetParent(node);
  is_first = false;
}
```

The `prev_old_fitch` and `prev_new_fitch` scratch buffers track the previous node's
original and modified Fitch sets across iterations, replacing the indexed lookups that
required knowing the previous node's identity.

### 3. Pre-allocated SrcRemovalResult

Added an in-place overload `ComputeInitialRemoval(src, SrcRemovalResult& out)` that writes
into an existing object instead of returning a new one. `UpwardTraversal` uses
`scratch.removal` throughout:

- **Level 1:** Directly populates `scratch.removal` with src's Fitch sets
- **Level 1→2 transition:** `ComputeInitialRemoval(src, scratch.removal)` overwrites in place
- **Level 2+:** `PropagateRemovalUpward(node, scratch.removal)` updates in place (unchanged)

The original `ComputeInitialRemoval(src)` returning by value is kept as a convenience
wrapper for test code.

### 4. Pointer-based Fitch access

Added `GetFitchSetPtr(NodeId)` to `TreeIndex`, returning `const uint8_t*` to the start of
a node's Fitch set row:

```cpp
const uint8_t* GetFitchSetPtr(NodeId node) const {
  return &fitch_sets_[node.value * num_variable_sites_];
}
```

Used throughout the hot path instead of per-site `GetFitchSet(node, si)` calls. This avoids
recomputing `node.value * num_variable_sites_` on every site iteration. Applied to:
- `ComputeMoveScoreCached` (src, dst, removal old/new, and per-path-node Fitch)
- `ComputeInitialRemoval` (src, src_parent, sibling)
- `PropagateRemovalUpward` (current_lca)
- `ComputeLowerBound` (src)
- `UpwardTraversal` (src, for level-1 removal setup)

### 5. Early-exit check

After processing each path node's sites, check whether the accumulated `total` can
possibly become negative given the remaining path nodes:

```cpp
nodes_remaining--;
if (total > 0 && nodes_remaining > 0 &&
    total > static_cast<int>(n_sites) * nodes_remaining) {
  return total;
}
```

Each remaining path node can save at most 1 per site (a node going from cost=1 to cost=0).
If `total` already exceeds `n_sites * nodes_remaining`, no amount of savings can make the
move profitable.

**Effectiveness:** Limited on seedtree because most moves cluster at the same LCA level
and have short dst paths. More impactful on trees with deeper structure where many
destinations are far from the LCA.

### 6. Removed SearchSubtreeDirectly

The `SearchSubtreeDirectly` method (which used the uncached `ComputeMoveScore`) was dead
code — not called from any hot path or test. Removed to reduce maintenance surface.

## Backward Compatibility

All original public API signatures are preserved as convenience wrappers:

```cpp
// These allocate internally and delegate to the scratch-buffer overloads
SrcRemovalResult ComputeInitialRemoval(NodeId src) const;
int ComputeMoveScoreCached(NodeId src, NodeId dst, NodeId lca,
                           const SrcRemovalResult& removal) const;
```

Test 120 (cached vs independent scoring) continues to construct `SrcRemovalResult` manually
and call `ComputeMoveScoreCached` without scratch buffers — verified working.

## Test Results

All 8 tests pass (113-120):

| Test | Result | Notes |
|------|--------|-------|
| 113: TreeIndex construction | PASS | |
| 114: Enumerator finds moves | PASS | 10 moves found |
| 115: Incremental vs full Fitch | PASS | 10 verified, 0 mismatches |
| 116: DAG grows | PASS | 11→42 nodes |
| 117: Valid DAG | PASS | parsimony=16 |
| 118: test_5_trees | PASS | local optimum, 0 moves |
| 119: Pruned vs exhaustive | PASS | 10=10, 0 missing, 0 extra |
| 120: Cached vs independent scoring | PASS | 48 checks, 0 mismatches |

## Remaining Bottlenecks

1. **Enumeration still dominates** — At radius 64, enumeration is 56s vs ~5s for apply +
   merge. The algorithmic complexity is unchanged; only the constant factor improved.

2. **No pruning** — The `allele_union` infrastructure from Phase 2 is still not used for
   pruning. A valid lower bound remains the main prerequisite.

3. **Apply phase** — At ~5s per radius, `InitHypotheticalTree` is the next target if
   enumeration is optimized further.

4. **Debug build overhead** — These benchmarks use a debug build (`build-deb-mock`). A
   release build with `KEEP_ASSERTS=on` would give more realistic absolute timings while
   preserving correctness checks.

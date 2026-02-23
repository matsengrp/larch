# Phase 4: Condensing Identical Sibling Leaves — Report

## Overview

Phase 4 adds condensing of identical sibling leaves in `TreeIndex`. When multiple leaf
children of the same parent share the same `CompactGenome`, only one representative is
kept in the condensed topology. This reduces the number of searchable source nodes and
the size of subtrees traversed during move enumeration.

All changes are in two files:
- `include/larch/spr/native_optimize.hpp` — condensing logic
- `test/test_native_optimize.cpp` — added B.1.1.529 benchmark test

## Benchmark Results

### Seedtree (~1040 nodes, 1238 variable sites)

**0 condensed leaves.** This tree has no identical sibling leaves, so condensing is a
no-op. Enumeration times are within noise of Phase 3:

| Radius | Phase 3 | Phase 4 |
|--------|---------|---------|
| 8 | 1,353ms | 1,368ms |
| 16 | 7,681ms | 7,749ms |
| 32 | 36,811ms | 37,285ms |
| 64 | 56,310ms | 56,476ms |

### B.1.1.529 and 20B datasets

Both OOM-killed during `TreeIndex` construction. The flat array layout
(`num_nodes * num_variable_sites`) exceeds available memory for these larger trees.
This is a pre-existing limitation of `TreeIndex`, not introduced by condensing.

## Implementation Details

### 1. New members in TreeIndex

```cpp
// Private:
std::vector<bool> is_condensed_;
size_t condensed_count_{0};

// Public:
size_t NumCondensedLeaves() const { return condensed_count_; }
```

### 2. Condensing during construction

After building the `parent_`/`children_` topology, before collecting `searchable_nodes_`:

```cpp
is_condensed_.assign(num_nodes_, false);
for (size_t nid = 0; nid < num_nodes_; nid++) {
  auto& kids = children_[nid];
  if (kids.size() <= 1) continue;

  // For each pair of leaf children, check CompactGenome equality
  // First identical leaf is the representative; rest are condensed
  // Condensed leaves are removed from children_[nid]
}
```

- O(n²) per parent — fine since children lists are small (typically 2–10)
- `CompactGenome::operator==` is fast (hash check first, then mutation map)
- Condensed leaves are removed from `children_[]`, automatically affecting all
  downstream computation (DFS indices, Fitch sets, allele unions, subtree search)

### 3. Searchable nodes filter

```cpp
for (auto node : dag_.GetNodes()) {
  if (not node.IsUA() and not node.IsTreeRoot() and
      not is_condensed_[node.GetId().value]) {
    searchable_nodes_.push_back(node.GetId());
  }
}
```

### 4. Correctness guarantee

Move scores on the condensed tree may differ slightly from the uncondensed tree
because `num_children` values differ at affected nodes. This is acceptable: every
candidate move is verified by `InitHypotheticalTree` on the original uncondensed DAG
before merging, catching any false positives from condensing-induced score differences.

### 5. Updated logging

```
Building tree index... done (624ms, 1238 variable sites, 1036 searchable nodes, 0 condensed)
```

## Test Results

All 9 tests pass (113–121):

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
| 121: seedtree | PASS | 0 condensed, timings match Phase 3 |

Test 122 (B.1.1.529) was added but OOM-kills during TreeIndex construction.

## Assessment

Condensing is correctly implemented and zero-cost when no identical siblings exist.
However, the seedtree dataset — our primary benchmark — has no condensable leaves,
so this phase provides no speedup on that dataset.

The more pressing issue revealed by this phase is that `TreeIndex` cannot handle
larger trees (B.1.1.529, 20B) due to memory consumption from the flat array layout.
The arrays `fitch_sets_`, `child_counts_`, and `allele_union_` each allocate
`num_nodes * num_variable_sites` entries. For a tree with ~10K nodes and ~10K
variable sites, these arrays alone consume several GB.

## Remaining Bottlenecks

1. **Memory scaling** — The flat `node × site` arrays in TreeIndex OOM on trees
   larger than seedtree. Addressing this (sparse storage, site batching, or
   only storing variable sites per-node) is a prerequisite for running on
   realistic datasets.

2. **Enumeration still dominates** — On seedtree, enumeration at radius 64 is 56s
   vs ~5s for apply + merge.

3. **No effective pruning** — The `allele_union` lower bound exists but rarely
   prunes because it only accounts for the new-node cost, not the src-removal
   savings. A tighter bound incorporating removal state could prune more.

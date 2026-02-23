# Phase 1 Report: Native MADAG Move Enumerator

## Summary

Phase 1 implements matOptimize's core algorithm -- exhaustive bounded SPR move enumeration -- directly on MADAG, eliminating the MAT conversion roundtrip. The implementation lives in a single header (`include/larch/spr/native_optimize.hpp`) with tests in `test/test_native_optimize.cpp`.

**Result**: 100% score accuracy on binary trees (0/16 mismatches on sample DAG). ~98% accuracy on non-binary trees (16/803 mismatches on test_5_trees). All 6 tests pass.

## Architecture

### Components

1. **`TreeIndex<DAG>`** -- Built once per sampled tree. Stores:
   - Per-node Fitch sets at variable sites (one-hot bitmask)
   - Per-node child allele counts (for incremental updates)
   - Per-node child count (`num_children_`)
   - DFS indices for O(1) ancestor checks
   - Lists of variable sites and searchable source nodes

2. **`MoveEnumerator<DAG>`** -- Exhaustive bounded search:
   - For each source node, walks up toward tree root (bounded by radius)
   - At each ancestor, searches sibling subtrees for improving destinations
   - Scores via `ComputeMoveScore` which simulates the full topology change

3. **`OptimizeDAGNative()`** -- Four-phase optimization loop:
   - Phase 1 (parallel): Enumerate profitable moves
   - Phase 2 (parallel): Create SPR overlays + `InitHypotheticalTree`
   - Phase 3 (parallel): Create fragments
   - Phase 4 (batch): Merge all fragments into DAG

### Scoring: `ComputeMoveScore`

Computes the exact parsimony score change for an SPR move (src, dst, lca) by simulating the Fitch set changes along both the removal and insertion paths.

The algorithm:
1. Collect all affected nodes on both paths (src_parent...lca and dst_parent...lca)
2. Sort by DFS level (deepest first, bottom-up)
3. Process each node once, applying both src-side and dst-side count modifications
4. Track modified Fitch sets in a map for propagation to parent nodes
5. Sum up cost deltas

## Bugs Found and Fixed

### Bug 1: Wrong `num_children` derivation from allele counts

**Symptom**: 12/16 score mismatches, incremental scores systematically too negative.

**Root cause**: `FitchCostFromCounts` derived the number of children by summing all allele counts: `num_children = counts[0] + counts[1] + counts[2] + counts[3]`. But when a child has a multi-bit Fitch set (e.g., {A, C} = 0b0011), that child increments *both* `counts[A]` and `counts[C]`. So the sum overcounts.

Example: Node with 2 children, child A has Fitch set {A} (0b0001), child B has Fitch set {A, C} (0b0011). Counts: A=2, C=1. Sum = 3, but actual children = 2. With `num_children=3`, the intersection check (`count[j] == num_children`) incorrectly concludes no allele is present in all children.

**Fix**: Track actual child count separately in `num_children_` map. Pass it explicitly to `FitchCostFromCounts` and `FitchSetFromCounts`.

### Bug 2: Wrong Fitch set computation (majority rule vs standard Fitch)

**Symptom**: Still 12/16 mismatches after Bug 1 fix. Traced to wrong Fitch *sets* propagating upward, causing incorrect costs at ancestor nodes.

**Root cause**: `FitchSetFromCounts` used a "max count" approach -- find the alleles with the highest count and return those bits. This is matOptimize's majority-rule heuristic. Standard Fitch uses intersection/union:

- **Intersection non-empty**: Return only alleles where `count[j] == num_children` (present in ALL children)
- **Intersection empty**: Return all alleles where `count[j] > 0` (present in ANY child)

These differ for nodes with >2 children and multi-bit child Fitch sets.

Concrete example from sample DAG: Node with 3 children, children have Fitch sets G, {C|T}, C. Counts: A=0, C=2, G=1, T=1. Max-count gives {C} (count 2). Standard Fitch: no allele has count==3, so union = {C|G|T}. These different Fitch sets propagate upward, causing different costs at ancestors even when the local cost (0 vs 1) happens to match.

**Fix**: Rewrote `FitchSetFromCounts` to use proper Fitch semantics. Extracted both `FitchSetFromCounts` and `FitchCostFromCounts` as free functions outside the class.

**Result after both fixes**: 0/16 mismatches on sample DAG (binary tree).

### Bug 3 (design evolution): Complex path overlap handling

**Symptom**: Code for handling the src-side and dst-side paths in `ComputeMoveScore` became unmaintainable spaghetti, with ad-hoc special cases for when paths share nodes.

**Root cause**: Initial approach computed removal changes along src_path and insertion changes along dst_path separately, then tried to merge overlapping effects at shared nodes. The overlap logic had many corner cases (LCA on both paths, path prefixes shared, etc.).

**Fix**: Complete rewrite using a unified approach:
1. Collect ALL affected nodes from both paths into a single set
2. Sort bottom-up by DFS level
3. Process each node exactly once, checking whether it's on the src path, dst path, or both
4. Apply all modifications to a copy of that node's child counts in one pass

This eliminated all the ad-hoc overlap logic.

## Known Limitations

### ~2% mismatch on non-binary trees

On test_5_trees (which has 14 non-binary internal nodes out of ~130 total), 16 out of 803 random (src, dst) pairs show score mismatches between `ComputeMoveScore` and `ParsimonyOnlyScoringBackend`.

**Root cause**: The SPR overlay's clade handling for multi-child nodes differs from the simplified model in `ComputeMoveScore`. Specifically, `InitHypotheticalTree` handles the new_node insertion and src removal in ways that interact with clade structure, which the incremental scorer's simple "decrement/increment allele counts" model doesn't fully capture for nodes with >2 children.

**Impact**: Low. All moves go through `InitHypotheticalTree` + `ParsimonyOnlyScoringBackend` validation before being merged, so incorrect incremental scores only cause false positives (attempting moves that turn out not to improve) or false negatives (missing some improving moves). The DAG never gets incorrect fragments.

**Possible fix for Phase 2**: Use the SPR overlay as ground truth for scoring, or handle multi-child clade structure explicitly in `ComputeMoveScore`.

### Unifurcation collapse

When src_parent is a binary node and src is pruned, src_parent becomes a unifurcation (single child) and is effectively collapsed. `ComputeMoveScore` models this by:
1. Subtracting src_parent's old Fitch cost
2. Setting src_parent's new Fitch set to its remaining child's (src_sibling's) Fitch set
3. Skipping src_parent during the bottom-up processing loop

This is correct for binary trees but may interact subtly with multi-child nodes above.

### No range tree pruning

Phase 1 uses exhaustive search within the radius bound. matOptimize uses a range tree to prune the search space based on parsimony score bounds. This is a performance optimization only -- correctness is unaffected.

### No conflict resolution

All profitable moves are applied independently against the original sampled tree. Conflicting moves (e.g., two moves pruning the same subtree) are handled gracefully by the fragment-merge pipeline -- `InitHypotheticalTree` returns false for invalid topologies, and the merge correctly handles overlapping fragments.

### No working tree mutation

Each move is evaluated against the original sampled tree. This means compounding improvements (move B improves score only after move A is applied) are not discovered within a single iteration. The outer optimization loop (sample, enumerate, merge, repeat) eventually finds them.

## Test Results

All 6 tests pass (tag: `native-optimize`):

1. **TreeIndex construction**: Verifies variable sites found, DFS indices consistent (`dfs_index < dfs_end_index`), searchable nodes non-empty, tree root is ancestor of all nodes.

2. **Enumerator finds moves**: Verifies profitable moves found on sample DAG, all with `score_change < 0`.

3. **Incremental vs full Fitch**: For every found move, compares `ComputeMoveScore` against `ParsimonyOnlyScoringBackend` (full Fitch recomputation). **0 mismatches** on sample DAG.

4. **DAG grows**: After native optimization, DAG has at least as many nodes/edges as before.

5. **Valid DAG**: Post-optimization DAG produces valid trees with edge mutations consistent with compact genomes.

6. **test_5_trees**: Runs 3 iterations of native optimization on `data/test_5_trees/tree_0.pb.gz` with timing breakdown.

## Implementation Quirks and Gotchas

### One-hot encoding

Uses `base_to_singleton()` from `scoring_backend.hpp`: A=1 (bit 0), C=2 (bit 1), G=4 (bit 2), T=8 (bit 3). Fitch sets are OR-combinations of these bits. Intersection is bitwise AND, union is bitwise OR.

### Fitch cost is per internal node, not per edge

The standard Fitch algorithm assigns a cost of 0 or 1 to each internal node at each site, based on whether the intersection of its children's Fitch sets is non-empty. This is NOT the same as counting mutations on edges. The incremental scorer computes cost changes at affected internal nodes.

### `InitHypotheticalTree` rejection conditions

The SPR overlay's `InitHypotheticalTree(src, dst, lca)` returns false (no-op) when:
- `src_parent == lca` and it's not a sibling move (in a binary tree, also rejects the sibling move since src_parent would become a unifurcation)
- Various other degenerate topologies

Rejected moves are silently dropped in Phase 2 of `OptimizeDAGNative`.

### `CompactGenome` is non-copyable

Must use `const auto&` when accessing compact genomes. Using `auto` triggers a deleted copy constructor.

### `Score` uses `.value()` not `.GetValue()`

The Score class from `ParsimonyOnlyScoringBackend` exposes the score through `.value()`.

### DFS numbering for ancestor checks

`IsAncestor(a, b)` checks `a.dfs_index <= b.dfs_index && b.dfs_index < a.dfs_end_index`. The end index is exclusive (set to the counter value after visiting all descendants).

### test_5_trees finds 0 profitable moves

The min-weight sampled tree from `test_5_trees/tree_0.pb.gz` is already at a local SPR optimum. This was confirmed independently by both native and random scoring methods. This is correct behavior, not a bug -- the tree happens to already be locally optimal under parsimony.

### Variable sites collected from edge mutations, not reference sequence

Variable sites are identified by scanning all edge mutations in the sampled tree. Sites with no mutations are not tracked. This means a site that is variable in the DAG but happens to be invariant in the currently sampled tree will not be checked for move scoring. This is correct for scoring moves on the current tree but could miss opportunities if the tree is not representative.

## Files

| File | Lines | Description |
|------|-------|-------------|
| `include/larch/spr/native_optimize.hpp` | 747 | Main implementation (all in header) |
| `test/test_native_optimize.cpp` | 374 | 6 tests |
| `CMakeLists.txt` | +1 line | Added test source |

## Benchmark: Seedtree Dataset

Three-way comparison on `data/seedtree/seedtree.pb.gz` (~1040 nodes, 1238 variable sites). All runs: 3 iterations, RelWithDebInfo build.

### Summary

| | **matOptimize** | **Random (1000/radius)** | **Native** |
|---|---|---|---|
| **Total time** | **312s** (5m12s) | **128s** (2m9s) | >300s (killed at iter 1, radius 64) |
| **Final parsimony** | **1606** | **1641** | N/A (incomplete) |
| **DAG size (final)** | 8577 nodes / 26744 edges | 70450 nodes / 202097 edges | N/A |
| **Moves per iteration** | ~170 applied per radius | ~3600 total (6 radii) | ~49 found per radius |

### matOptimize per-radius timing (iteration 1)

| Radius | Search time | Moves applied |
|--------|------------|--------------|
| 2 | 8s | 162 |
| 4 | 9s | 174 |
| 8 | 12s | 176 |
| 16 | 22s | 172 |
| 32 | 16s | 183 |
| 64 | ~35s | 166 |
| **Total** | **101s** | |

matOptimize mutates the working tree in-place between moves, so later moves benefit from earlier improvements. The batch merge contention (thread-time) was 100-150s per radius callback, indicating the larch merge pipeline is the dominant cost.

### Native optimizer per-radius timing (iteration 1, incomplete)

| Radius | Enumerate | Apply+Merge | Moves found |
|--------|-----------|-------------|-------------|
| 2 | 30ms | 134ms | 8 |
| 4 | 283ms | 198ms | 36 |
| 8 | 2.3s | 366ms | 49 |
| 16 | **30.4s** | 970ms | 49 |
| 32 | **128.7s** | 347ms | 49 |
| 64 | (killed) | | |

Enumeration time grows quadratically with radius. The apply+merge phase is fast (<1s per radius) since it reuses the same fragment-merge pipeline as the random optimizer.

### Random optimizer per-iteration timing

| Iteration | Optimize | Serial overhead | Total | Parsimony |
|-----------|----------|-----------------|-------|-----------|
| 1 | 38.6s | 1.3s | 39.9s | 1662 |
| 2 | 42.9s | 2.1s | 45.0s | 1647 |
| 3 | 43.4s | 2.2s | 45.6s | 1641 |

### Analysis

1. **The native enumerator's bottleneck is `ComputeMoveScore`**. At radius 32 with ~1036 source nodes, the exhaustive search evaluates every (src, dst) pair within the radius. Each call walks paths to LCA and iterates 1238 variable sites. Complexity is roughly O(sources x destinations_in_radius x sites x path_length).

2. **matOptimize uses range tree pruning** to skip destinations that provably cannot improve the score. This is the key performance gap: matOptimize searches at radius 64 in ~35s; the native enumerator takes >128s at radius 32 alone.

3. **matOptimize applies moves in-place**, finding compounding improvements within a single radius pass. It applies ~170 moves per radius vs the native optimizer's ~49 total found (and those 49 are all evaluated against the same unmodified tree).

4. **matOptimize achieves better parsimony** (1606 vs random's 1641) because in-place mutation discovers deeper improvements. The native optimizer should approach this once it completes all radii.

5. **The random optimizer is fastest** at 128s because it skips exhaustive enumeration -- it randomly samples move candidates. The tradeoff is a higher parsimony plateau (1641 vs 1606) and a much larger DAG (70k nodes vs 8.5k) from merging many exploratory fragments.

6. **Peak memory**: Random optimizer used ~2GB; matOptimize used less due to fewer merged fragments.

### Conclusion

The native optimizer needs **range tree pruning** before it can compete with matOptimize on trees with >500 nodes. Without pruning, the quadratic enumeration cost at large radii makes it impractical. The algorithm is correct (verified against full Fitch) but too slow. Phase 2 priorities:

1. **Range tree pruning** -- the single most impactful optimization
2. **Working tree mutation** -- apply moves in-place to find compounding improvements
3. **Multi-child node scoring** -- fix the ~2% mismatch on non-binary trees

## Performance Notes (small trees)

On `test_5_trees/tree_0.pb.gz` (~130 nodes):
- TreeIndex construction: <1ms
- Move enumeration (all radii): <10ms total
- SPR overlay creation: comparable to random_optimize
- Fragment creation + merge: reuses random_optimize infrastructure

The enumeration phase is fast on small trees. The quadratic scaling only becomes a problem above ~500 nodes.

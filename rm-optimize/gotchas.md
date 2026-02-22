# Gotchas, Edge Cases, and Non-Obvious Behaviors

## matOptimize Gotchas

### 1. Radius is packed with flags

The `radius` parameter is not just a number. It's a bitmask:
```cpp
#define RADIUS_MASK  0x3fffffff
#define DRIFT_MASK   0x80000000
#define ALL_DIR_MASK 0x40000000
```
Extract actual radius with `radius & RADIUS_MASK`. Check drift mode with
`radius & DRIFT_MASK`. Easy to forget when porting.

### 2. Mutations_Collection is sorted by position

All mutation vectors are sorted by genomic position. This is not just a
convention -- the merge algorithms use two-pointer techniques that assume sorted
order. Inserting unsorted mutations will silently corrupt state.

MADAG equivalent: `ContiguousMap` is also sorted. This matches well.

### 3. Negative positions mean masked mutations

In MAT::Mutation, a negative `position` indicates a masked (ignored) mutation.
These are skipped during scoring but preserved during moves. Don't filter them
out prematurely.

### 4. par_mut_nuc packing

Parent and mutant nucleotides are packed into a single byte:
```cpp
par_nuc = (par_mut_nuc >> 4) & 0xf
mut_nuc = par_mut_nuc & 0xf
```
Easy to get the shift direction wrong. Also, one-hot encoding means A=1, C=2,
G=4, T=8 -- NOT 0,1,2,3.

### 5. boundary1_all_major_allele packing

Same byte packs two different sets:
```cpp
all_major_allele = boundary1_all_major_allele & 0xf     // lower nibble
boundary1_allele = (boundary1_all_major_allele >> 4) & 0xf  // upper nibble
```

### 6. The "sensitive increment" is not just boundary1

An allele is sensitive if incrementing it could change the major allele set.
This includes:
- Alleles already in the major set (incrementing preserves majority)
- Alleles in boundary1 (incrementing ties with major)
- But NOT alleles below boundary1 (incrementing still doesn't reach major)

### 7. move_node() uses INVERT_MERGE for downward mutations

When merging mutations from dst toward dst's parent, the merge direction is
inverted. This is because we're conceptually moving UP the tree (from dst's
perspective), so parent and child roles swap. Using normal merge here will
produce wrong mutations.

### 8. Delayed deletion is critical

After `move_node()`, nodes marked for deletion are NOT immediately freed.
Other moves in the same batch may still reference them through
`potential_crosses` or the conflict resolver. Deleting early causes
use-after-free.

### 9. reassign_states is expensive and sometimes unnecessary

`reassign_states()` does a full tree Fitch-Sankoff reconstruction. It's called
after large batch operations but NOT after every individual move. Over-calling
it wastes time. Under-calling it can leave inconsistent state.

### 10. Condensing changes node count

After condensing, the tree has fewer nodes. All node indices, DFS ordering,
and range trees become invalid. Must rebuild after condensing. Similarly,
uncondensing invalidates everything.

## MADAG / Larch Gotchas

### 11. CompactGenome is relative to reference, not parent

`CompactGenome` stores mutations from the **Universal Ancestor's reference
sequence**, not from the parent node. To get parent→child mutations, use
`ToEdgeMutations(parent_cg)` or compare the two compact genomes.

This is different from MAT where `node->mutations` are parent-relative.

### 12. EdgeMutations entries have (parent_base, child_base) pairs

Each entry in `EdgeMutations` is a `std::pair<MutationBase, MutationBase>` at a
`MutationPosition`. The first element is the parent's base, the second is the
child's base. Getting these swapped silently produces wrong results.

### 13. DAG vs Tree: parent count matters

In a DAG, nodes can have multiple parents. In a sampled tree, each node has
exactly one parent (except UA which has none). Code that assumes single-parent
will break on DAGs. Always sample a tree first for optimization.

### 14. BuildConnections() must be called

After constructing or modifying a MADAG, `BuildConnections()` must be called
to establish parent-child relationships from the stored edge data. Forgetting
this leaves the DAG in an unusable state with no traversal possible.

### 15. ComputeCompactGenomes() is order-dependent

Must be called after `BuildConnections()` and after edge mutations are set.
Calling it before edges are populated produces zero-mutation compact genomes.

### 16. Fragment edge IDs are not contiguous

`MakeFragment()` produces sparse node and edge IDs. They correspond to IDs in
the original tree, not sequential 0..N. Code that assumes contiguous IDs will
index out of bounds.

### 17. CollapseEmptyFragmentEdges can remove "important" nodes

If a node has the same compact genome as its parent and isn't a leaf or
move-related node, it gets collapsed. This is correct for merging but means
the fragment doesn't represent the full subtree topology.

### 18. MappedNodes must be overlaid before merge

Fragments need `Deduplicate<MappedNodes>` overlaid before being passed to
`merge.AddDAG()`. Without this, the merge can't track which nodes in the
fragment correspond to existing DAG nodes.

### 19. The UA node is special

The Universal Ancestor (root of all roots) has no parent and no compact genome
mutations. Its compact genome is empty (all positions equal to reference).
SPR moves must never use UA as src or dst. The random move generator explicitly
checks for this.

### 20. IsDescendant walks parent pointers

`RandomMoveGenerator::IsDescendant()` walks from descendant toward root
checking for ancestor. This is O(depth), not O(1). For frequent ancestor
checks during enumeration, precompute DFS indices instead.

## Performance Gotchas

### 21. Range tree memory

One range tree per variable site. For a genome with 30k sites and 100k nodes,
this can be significant. matOptimize uses FIFO allocators to manage this.
Consider memory pools or arena allocation for the MADAG version.

### 22. ContiguousMap merge is O(n+m)

Merging two ContiguousMaps (sorted vectors) is O(n+m) with a two-pointer scan.
This is efficient for similarly-sized maps but wasteful when inserting a single
element into a large map. Consider keeping a separate "pending insertions"
buffer and batch-merging periodically.

### 23. Hash map allocation in Merge

`Merge::result_nodes_` and `result_edges_` use `GrowableHashMap`. For large
DAGs with millions of nodes, hash map operations become a bottleneck. The
batched merge (`AddDAGs`) amortizes this.

### 24. TBB vs Taskflow

Larch is migrating from TBB to Taskflow. matOptimize uses TBB's flow graph
extensively. The native implementation should use Taskflow (or at minimum
larch's `ParallelForEach`) to stay consistent with the migration direction.

### 25. DISABLE_PARALLELISM affects TBB globals

When built with `DISABLE_PARALLELISM=yes`, `tbb::global_control` limits
parallelism to 1. This affects ALL TBB usage in the process, including
matOptimize's internal TBB code. If the native implementation uses Taskflow,
it needs its own parallelism control separate from TBB.

## Correctness Gotchas

### 26. Score change sign convention

Negative score_change = improvement (lower parsimony). Zero = neutral (accepted
in drift mode). Positive = worse (rejected). Be consistent.

matOptimize uses this convention. Larch's `ParsimonyOnlyScoringBackend` also
uses `new_score - old_score` (negative = better).

### 27. Move validity after tree mutation

After applying a move to the working tree, previously-found moves may become
invalid:
- src or dst may have been absorbed (unifurcation cleanup)
- The path between src and dst may have changed
- Score predictions may be wrong

Always re-validate deferred moves on the current tree state.

### 28. Fitch set vs major allele: subtle difference

The Fitch set at a node is the set of alleles that achieve optimal parsimony
for the subtree. The major allele is the most frequent allele among children.
These are NOT the same:

- Fitch set considers the full subtree (all descendants)
- Major allele only considers direct children
- For binary trees with no polytomies, they often agree
- For polytomies, they can differ

matOptimize tracks major alleles (simpler, incremental). Larch's backend
computes Fitch sets (correct, but expensive). For the native implementation,
major allele tracking is sufficient for incremental scoring.

### 29. Empty edge mutations after collapse

`CollapseEmptyFragmentEdges()` can produce fragments where some expected edges
are missing. The merge system handles this correctly, but debugging tools that
expect complete subtree topology will be confused.

### 30. Multiple trees in the DAG may share the same parsimony score

When the DAG grows, it accumulates multiple trees with the same optimal
parsimony score. `MinWeightSampleTree()` picks one arbitrarily. Different
samples may lead to different optimization paths. This is expected behavior,
not a bug.

# Conflict Resolution and Move Scheduling

## Why Conflicts Arise

SPR moves can be found in parallel across different source nodes. But two moves
may affect overlapping regions of the tree -- applying both would invalidate
the score computation of at least one.

Two moves **conflict** when their "paths" intersect. A move's path includes all
nodes from src to LCA and from dst to LCA, plus all ancestors whose major
allele changes.

## matOptimize's Conflict Resolution

### Data Structure: potential_crosses

An array indexed by BFS index. Each entry points to the "best move so far" that
crosses through that node.

```cpp
// Conceptually:
potential_crosses[bfs_idx] = {move_ptr, score_change}
```

### Resolution Algorithm

For each new candidate move:

1. **Collect conflicts**: Walk the move's path. At each node, check
   `potential_crosses[node.bfs_idx]`. If non-null, record the conflicting move.

2. **Find minimum competitor**: Among all conflicting moves, find the one with
   the best (most negative) score_change.

3. **Win/lose decision**:
   - If new move's score < best competitor's score: new move wins
   - Otherwise: new move is deferred

4. **If new move wins**:
   - Remove all losing moves from `potential_crosses`
   - Register new move at all nodes on its path
   - Losing moves go to the `deferred` list

5. **If new move loses**: Add to `deferred` list immediately

### Why BFS ordering?

BFS indices ensure parents have lower indices than children. This provides a
consistent ordering for conflict resolution without race conditions.

## Deferred Move Recycling

After the initial search-and-resolve phase, deferred moves get a second chance:

```
while deferred_moves is not empty:
    clear potential_crosses
    for each deferred move (src, dst):
        re-evaluate on the MODIFIED tree
        if still profitable and non-conflicting:
            accept
        elif conflicts with better move:
            defer again
    apply accepted moves
    repeat
```

This recycling continues until no more deferred moves can be accepted or the
improvement plateaus.

### Why recycling matters

Consider three moves A, B, C where A conflicts with both B and C:
- Initial resolution: A wins (best score), B and C deferred
- After applying A, the tree changes
- B's path may no longer conflict with anything → B can now apply
- C's path may still conflict with B → C deferred again

Recycling extracts maximum value from each search pass.

## Move Scheduling

Non-conflicting moves are ordered for application. The order matters because
`move_node()` modifies the tree, and subsequent moves need valid node pointers.

matOptimize's approach:
- Process moves in reverse DFS order (deepest source nodes first)
- This ensures parent nodes are still valid when their children are moved
- Delayed deletion prevents use-after-free

## Parallel Safety

The conflict resolution itself runs **sequentially** (single resolver node in
the TBB flow graph). This is intentional:
- Conflict checking requires a consistent view of `potential_crosses`
- Making it parallel would require locks at every node

The search phase runs in parallel (multiple searcher nodes). Results flow into
the single resolver.

## Mapping to MADAG

### What larch currently does

Larch's random move optimization has no conflict resolution because moves are
evaluated independently against the original tree. Each move produces an
independent fragment that gets merged. Conflicting fragments are fine -- the
merge handles them by creating alternative paths in the DAG.

### Options for native implementation

#### Option A: No conflict resolution (DAG-native approach)

Evaluate all moves independently. Merge all improving fragments into the DAG.
The DAG naturally represents all alternatives. Re-sample and repeat.

**Pros**: Simple, naturally parallel, leverages DAG structure
**Cons**: Doesn't find compounding improvements, may produce large DAGs with
many suboptimal paths

#### Option B: matOptimize-style conflict resolution

Maintain a working tree. Find moves in parallel. Resolve conflicts. Apply to
working tree. Recycle deferred moves.

**Pros**: Finds better individual trees, compounding improvements
**Cons**: More complex, requires working tree mutation

#### Option C: Batch-and-merge hybrid

1. Find all improving moves in parallel (no conflict resolution)
2. Group non-conflicting moves (using path intersection check)
3. Apply each non-conflicting group as a batch
4. For each batch, create fragments and merge into DAG

This gets some compounding benefit without full conflict resolution.

### Conflict detection on MADAG

For MADAG, conflict detection needs:
- A way to identify a move's "path" (nodes affected by the move)
- A way to check if two paths intersect

Using NodeId instead of BFS index:
- Could use a hash set of affected NodeIds per move
- Two moves conflict if their affected sets intersect
- More memory but simpler than maintaining a BFS-indexed array

Alternative: Compute DFS indices on the sampled tree (fast, O(N)) and use the
same BFS-indexed approach as matOptimize.

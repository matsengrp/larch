# matOptimize Algorithm Overview

## The Big Picture

matOptimize iteratively improves a phylogenetic tree's parsimony score by
finding and applying SPR (Subtree Prune and Regraft) moves. An SPR move takes a
subtree rooted at a **source** node, detaches it, and reattaches it as a sibling
of a **destination** node. If the resulting tree has a lower parsimony score
(fewer total mutations across all edges), the move is an improvement.

## Three-Level Loop Structure

### Level 1: Radius Iteration (main.cpp)

The outermost loop controls the **search radius** -- the maximum distance (in
tree edges) between source and destination nodes.

```
for radius in [initial, initial*2, initial*4, ...]:
    nodes_to_search = find_nodes_to_move(tree)
    optimize_inner_loop(nodes_to_search, tree, radius, callback)
    if no_improvement and radius < max:
        double radius, reset search
    if no_improvement at max radius:
        enter drift mode or terminate
```

**Radius evolution:**
- Starts small (e.g., 2) and doubles after each pass with no improvement
- Capped at `2 * tree.max_level` (diameter of the tree)
- Stored as a bitmask with flags: `RADIUS_MASK | DRIFT_MASK | ALL_DIR_MASK`

**Why radius matters:** Small radii find local improvements cheaply. Large radii
find global rearrangements but are expensive. Doubling explores local first,
then expands.

### Level 2: Inner Loop (optimize_inner_loop.cpp)

Within a fixed radius, the inner loop processes batches of source nodes:

```
while nodes_to_search is not empty:
    shuffle(nodes_to_search)       # randomize processing order
    optimize_tree_main_thread(nodes_to_search, tree, radius, callback)
    if parsimony_score unchanged:
        break                      # no improvement at this radius
    # otherwise loop to find more moves with updated tree
```

Each call to `optimize_tree_main_thread` may find and apply multiple moves,
changing the tree. The loop continues until a full pass yields no improvement.

### Level 3: Search + Apply (optimize_tree.cpp)

The core search-and-apply pipeline uses a TBB flow graph:

```
fetcher → searcher (N copies) → conflict_resolver
                                       ↓
                                 schedule_moves
                                       ↓
                                  apply_moves
                                       ↓
                                 recycle_deferred
```

1. **Fetch**: Distribute source nodes to worker threads
2. **Search**: For each source, call `find_moves_bounded()` to enumerate
   profitable SPR moves within the radius
3. **Resolve conflicts**: Moves whose paths overlap cannot both apply. Keep the
   best one, defer the rest.
4. **Schedule**: Order non-conflicting moves for application
5. **Apply**: Execute moves on the tree (see move-application.md)
6. **Recycle**: Re-evaluate deferred moves on the now-modified tree

## Termination Conditions

The algorithm stops when ANY of:
1. No parsimony improvement for a full radius sweep
2. Improvement ratio drops below `min_improvement` threshold
3. `max_round` iterations completed
4. Drift iterations exhausted (when allowing equal-score moves)
5. External interrupt (SIGUSR2)

## Drift Mode

After no-improvement termination, the algorithm optionally enters **drift mode**:
- Accepts moves with `score_change == 0` (same parsimony)
- Purpose: escape local optima by exploring the neutral landscape
- Controlled by `drift_iterations` parameter
- Uses `DRIFT_MASK` flag in the radius variable

## Node Selection: find_nodes_to_move()

Not every node is worth searching as a source. The selection considers:
- All non-root, non-leaf internal nodes are candidates
- After the first pass, nodes near recently-changed regions are prioritized
- The `ALL_DIR_MASK` flag resets this to search everything

## Integration with Larch

In larch's optimization loop (`test_spr.cpp` / `larch-usher.cpp`):

```
for each iteration:
    1. merge.ComputeResultEdgeMutations()
    2. Sample min-weight tree from DAG
    3. Convert MADAG → MAT::Tree (BuildMAT)
    4. Condense the MAT
    5. Call optimize_inner_loop() with callbacks
    6. Callbacks receive moves, create SPR overlays, make fragments
    7. Fragments merged into DAG via merge.AddDAGs()
    8. Convert optimized MAT → MADAG (BuildFromMAT)
    9. Merge optimized tree into DAG
```

The MAT↔MADAG conversion is what we want to eliminate. The algorithm itself
(search, score, apply, conflict-resolve) needs to work directly on MADAG.

## Key Insight: Two Levels of Move Application

matOptimize applies moves **to the MAT** (mutating tree structure). Larch also
creates SPR overlays and fragments **from the MADAG**. These are conceptually
different:

- **MAT application**: Destructive mutation of the tree for continued search
- **MADAG fragments**: Non-destructive overlays merged into the growing DAG

For native MADAG optimization, we need both: mutate a working tree (for
continued search within a radius) AND produce fragments (for DAG growth).

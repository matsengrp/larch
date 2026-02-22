# Move Application: Backward and Forward Passes

After profitable moves are selected and conflicts resolved, they must be
applied to the tree. This is where the tree structure and node states are
actually modified.

## Overview

Move application has three stages:
1. **move_node()**: Physically rearrange the tree topology
2. **Backward pass**: Recompute major alleles bottom-up
3. **Forward pass**: Propagate state changes top-down

## Stage 1: move_node() -- Topology Change

Located in `apply_move/move_node_no_reassign.cpp`.

### What happens

1. **Detach src** from its parent:
   - Remove src from parent's children list
   - If parent now has 0 or 1 child, mark for deletion/absorption

2. **Accumulate mutations** along the path from src to dst:
   - Merge mutations upward from src to src_parent (using `KEEP_SELF`)
   - Merge mutations downward from dst to dst_parent (using `INVERT_MERGE`)
   - This ensures src's absolute state (root → src path) is preserved

3. **Insert src at dst**:
   - Split the edge from dst to dst_parent
   - Create new internal node between dst and dst_parent
   - src and dst become children of the new node
   - Call `merge_new_node_mutations()` to distribute mutations correctly

### Mutation Distribution at the Split

When creating a new internal node N with children src and dst:

```
Before:                After:
  P                      P
  |                      |
  dst                    N
                        / \
                      src  dst
```

Mutations on edge P→dst must be redistributed:
- **Shared mutations** (same allele in both src and dst): go on edge P→N
- **dst-only mutations**: go on edge N→dst
- **src-only mutations**: go on edge N→src
- If src and dst disagree: new node N gets the binary Fitch optimal choice

### Node Cleanup

After moving, some nodes may become redundant:
- **Unifurcation**: Parent of src's old position now has 1 child → absorb
  (merge parent and child into one node, concatenate mutations)
- **Empty**: Parent has 0 children → delete entirely

These are tracked in an `altered_nodes` list for the backward pass.

### State Invariant

The critical invariant: **every leaf's absolute genome is unchanged.** The path
of mutations from root to each leaf must produce the same sequence before and
after the move. `move_node()` achieves this by carefully merging mutations.

## Stage 2: Backward Pass -- Major Allele Recomputation

Located in `apply_move/backward_pass.cpp`.

### Why it's needed

After topology changes, nodes have new children. Their major allele sets
(most frequent allele among children) are stale. The backward pass fixes this
bottom-up.

### Algorithm

Uses a **max-heap** ordered by DFS index (processes deeper nodes first):

```
heap = all altered_nodes (nodes whose children changed)
while heap is not empty:
    node = heap.pop()  // deepest first
    old_major = node.major_allele_at_each_position
    recompute major alleles from node's current children
    if major_allele changed at any position:
        record state change in Altered_Node_t
        heap.push(node.parent)  // propagate upward
```

### Major Allele Computation

For binary nodes (`one_level_fitch_sankoff_binary.cpp`):
```
intersection = child1_allele & child2_allele
if intersection:
    major = intersection  (both children agree)
else:
    major = child1_allele | child2_allele  (disagree, keep both)
```

For polytomy nodes (`one_level_fitch_sankoff.cpp`):
```
count allele frequencies across all children
major = alleles with maximum count
boundary1 = alleles with count = max_count - 1
```

### Propagation

Changes propagate upward until they stabilize. In the best case, a change at
a leaf only affects a few ancestors. In the worst case, it reaches the root.

### Output

A list of `Altered_Node_t` entries, each containing:
- Node pointer
- List of (position, old_state, new_state) triples
- These feed into the forward pass

## Stage 3: Forward Pass -- State Propagation

Located in `apply_move/forward_pass.cpp`.

### Why it's needed

The backward pass changes which allele a node "prefers" (major allele). But the
node's children may still be using the old preference. The forward pass updates
children to follow the new parent state when possible.

### Algorithm

Uses a **min-heap** ordered by DFS index (processes shallower nodes first):

```
heap = all nodes with state changes from backward pass
while heap is not empty:
    node = heap.pop()  // shallowest first
    for each child of node:
        for each position where node's state changed:
            if child can follow new parent state (allele in child's Fitch set):
                update child's allele to match parent
                record child's state change
        if child's state changed:
            heap.push(child.parent)  // propagate downward
```

### The "Lazy Following" Principle

A node prefers to follow its parent's allele when possible, because this avoids
a mutation on the connecting edge. The forward pass enforces this greedily.

### Interaction with Backward Pass

The backward and forward passes together implement a local Fitch-Sankoff
reconstruction limited to the changed region of the tree. This is much cheaper
than a full tree reconstruction.

## matOptimize's reassign_states()

Separate from the incremental backward/forward passes, `reassign_states.cpp`
does a full tree Fitch-Sankoff reconstruction:
- Complete bottom-up pass over all nodes
- Complete top-down pass over all nodes
- Used after large structural changes or for verification
- Much more expensive than incremental passes

## Mapping to MADAG

### What larch already has

Larch's `InitHypotheticalTree()` + `MakeFragment()` does something analogous
but non-destructively:

1. Creates an overlay over the sampled tree
2. Computes new CompactGenomes for affected nodes (via Fitch)
3. Extracts changed nodes/edges as a fragment
4. Fragment is merged into the DAG

This approach does NOT modify the sampled tree. Each SPR move produces an
independent fragment.

### The key difference

matOptimize **mutates** the tree between moves. This means:
- Move N sees the tree as modified by moves 1..N-1
- Moves can build on each other within a radius
- The tree improves incrementally

Larch's current approach evaluates each move against the **original** sampled
tree. Moves don't see each other's effects. This is correct for DAG growth but
misses compounding improvements.

### Options for native implementation

1. **Stateless (current larch approach)**: Each move evaluated independently on
   the sampled tree. Simple, naturally parallel, but misses compound effects.

2. **Stateful (matOptimize approach)**: Maintain a working copy of the tree.
   Apply accepted moves to the working copy. Evaluate subsequent moves on the
   modified tree. More complex but finds better solutions.

3. **Hybrid**: Evaluate moves statelessly for a batch, apply the best
   non-conflicting moves to a working copy, then re-evaluate remaining moves on
   the updated tree.

The stateful approach requires implementing the backward/forward pass machinery
over MADAG's data structures (`CompactGenome`, `EdgeMutations`). The stateless
approach can use the existing `ParsimonyOnlyScoringBackend` but needs the
bounded search for efficiency.

## Node State Flags (matOptimize)

matOptimize tracks which nodes need reprocessing:
```cpp
SELF_CHANGED      // This node's mutations changed
ANCESTOR_CHANGED  // Parent's state changed (may affect this node)
DESCENDENT_CHANGED // Child's state changed (may affect this node)
SELF_MOVED        // This node was physically moved
```

These flags guide the backward and forward passes to only visit nodes that
actually need updates. For a MADAG implementation, similar dirty-tracking would
be needed for the stateful approach.

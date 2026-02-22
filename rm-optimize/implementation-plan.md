# Implementation Plan: Native MADAG Optimization

## Goal

Replace the matOptimize dependency with a native implementation that operates
directly on MADAG, eliminating the MAT↔MADAG conversion roundtrip.

## Architecture Decision: What to Build

### Recommended: Phased approach

Build in layers, each independently useful:

1. **Efficient move enumeration** on MADAG (the big win)
2. **Conflict resolution** for parallel move application
3. **Working tree mutation** for compounding improvements
4. **Condensing** for performance on large datasets

Each layer can be tested independently using the existing random moves
infrastructure as a baseline.

## Phase 1: Efficient Move Enumeration

### 1.1 Per-Node State: Allele Counts

Add a new Feature to MADAG nodes that tracks allele frequency information at
each variable site:

```cpp
struct AlleleCount {
    MutationPosition position;
    uint8_t major_allele;     // One-hot bitmask of most frequent allele(s)
    uint8_t boundary1_allele; // One-hot bitmask of alleles at count = max-1
    uint8_t allele_counts[4]; // Counts for A, C, G, T
};

// Per node: vector<AlleleCount> at variable sites only
```

**Initialization**: Bottom-up pass over the sampled tree, counting alleles at
each site from children's compact genomes.

**MADAG equivalent of MAT::Mutation**: A node's `CompactGenome` gives its
absolute state. The `EdgeMutations` on parent edge give the parent→child
mutations. Together, these provide all information that MAT::Mutation stores.

### 1.2 DFS Index Precomputation

Add DFS indices to the sampled tree:

```cpp
struct DFSInfo {
    size_t dfs_index;      // Preorder position
    size_t dfs_end_index;  // End of subtree in preorder
    size_t level;           // Depth from root
};
```

Compute once per sampled tree with a simple DFS traversal. Store as a vector
indexed by NodeId.

This enables O(1) ancestor checks:
```cpp
bool IsAncestor(NodeId a, NodeId b) {
    return dfs[b].dfs_index >= dfs[a].dfs_index &&
           dfs[b].dfs_index <= dfs[a].dfs_end_index;
}
```

### 1.3 Range Tree Construction

Build the range tree from MADAG's edge mutations:

```cpp
struct RangeTreeNode {
    uint32_t dfs_start, dfs_end;
    uint8_t min_level[4];  // Min depth per allele
    // Children indices in flat array
};

// One range tree per variable site
std::vector<std::vector<RangeTreeNode>> range_trees;
```

Construction:
1. For each variable site, collect all nodes with a mutation at that site
   (from edge mutations)
2. Sort by DFS index
3. Build hierarchical intervals with min_level per allele

### 1.4 Incremental Score Computation

Implement the upward-downward traversal using MADAG data:

```cpp
struct MutationCountChange {
    MutationPosition position;
    uint8_t decremented_allele;  // One-hot
    uint8_t incremented_allele;  // One-hot
    uint8_t par_state;           // One-hot
};
```

**Upward phase**: For source node `src`, walk toward root:
```
At each ancestor A:
    For each mutation M on edge A→child_toward_src:
        Update allele count at A: decrement child's allele
        Check if major allele changed
        Compute removal score change
    Search A's other children (downward phase) for destinations
```

**Downward phase**: For each sibling subtree of the current LCA:
```
At each node N in DFS order:
    Compute split cost: inserting src as sibling of N
    Check range tree: can going deeper improve?
    Update lower bound
    If lower_bound > best_score: prune
```

### 1.5 Candidate Evaluation Functions

Three functions to implement:

```cpp
int ComputeRemovalScore(NodeId src, const AlleleCounts& counts);
int ComputeSplitScore(NodeId src, NodeId dst, const AlleleCounts& counts);
int ComputeAboveLCAScore(NodeId lca, const AlleleCounts& old_counts,
                          const AlleleCounts& new_counts);
```

These correspond to matOptimize's `check_move_profitable_LCA()` and
`check_move_profitable_dst_not_LCA()`.

### Testing Phase 1

- Use `ParsimonyOnlyScoringBackend` as ground truth
- For each move found by the new enumerator, verify score with full Fitch
- Compare with random moves: new enumerator should find all improving moves
  that random sampling finds, plus more

## Phase 2: Conflict Resolution

### 2.1 Move Path Tracking

For each candidate move, record all affected NodeIds:

```cpp
struct MoveCandidate {
    NodeId src, dst, lca;
    int score_change;
    std::vector<NodeId> affected_nodes;  // Path from src to lca to dst
};
```

### 2.2 Conflict Detection

Two moves conflict if their `affected_nodes` sets intersect.

Options:
- **Hash set intersection**: O(path_length) per pair
- **BFS-indexed array** (matOptimize approach): O(path_length) total, amortized
- **Bitset per move**: O(N/64) per pair, fast with SIMD

For moderate tree sizes (up to ~100k nodes), the BFS-indexed array approach
scales well. Allocate once per sampled tree.

### 2.3 Priority Resolution

Same algorithm as matOptimize:
1. Check all nodes on the new move's path
2. If any node has a registered move with better score, defer new move
3. If new move is better, deregister all conflicting moves, register new move
4. Apply all registered (non-conflicting) moves

### 2.4 Move Recycling

After applying a batch, re-evaluate deferred moves on the modified tree.
This requires Phase 3 (working tree mutation) to be effective.

## Phase 3: Working Tree Mutation

### 3.1 Mutable Sampled Tree

Instead of using the sampled tree as read-only (current approach), maintain a
mutable working copy:

```cpp
auto working_tree = weight.MinWeightSampleTree({});
// Apply accepted moves to working_tree
// Re-evaluate remaining moves on working_tree
```

### 3.2 Backward Pass on MADAG

After applying a move to the working tree, update allele counts bottom-up:

```cpp
void BackwardPass(WorkingTree& tree, const std::vector<NodeId>& altered_nodes) {
    // Process in leaf-to-root order (reverse DFS)
    for each altered_node in reverse DFS order:
        recompute allele counts from children
        if allele counts changed:
            mark parent as altered
}
```

Allele counts are stored in the per-node `AlleleCount` structures from Phase 1.

### 3.3 Forward Pass on MADAG

Propagate state changes top-down:

```cpp
void ForwardPass(WorkingTree& tree, const std::vector<NodeId>& changed_nodes) {
    // Process in root-to-leaf order (DFS order)
    for each changed_node in DFS order:
        for each child:
            if child can follow new parent state:
                update child's compact genome
                mark child as changed
}
```

### 3.4 Fragment Extraction

After mutating the working tree, the existing `MakeFragment()` infrastructure
can extract changed regions for merging into the DAG. But we need to track
which nodes changed during working tree mutation.

Alternative: After all moves in a radius, compare working tree to original
sampled tree and produce a single large fragment (or the whole tree).

## Phase 4: Condensing

### 4.1 Leaf Condensing on MADAG

Group leaf siblings with identical compact genomes:

```cpp
void CondenseLeaves(WorkingTree& tree) {
    for each internal node:
        group = children with same CompactGenome that are leaves
        if group.size() > 1:
            keep one representative, mark others as condensed
            store mapping: representative → [condensed leaves]
}
```

### 4.2 Uncondensing

After optimization, restore condensed leaves:

```cpp
void UncondenenseLeaves(WorkingTree& tree) {
    for each condensed group:
        restore all leaves as siblings of the representative
}
```

## Concept Mapping: matOptimize → MADAG

| matOptimize Concept | MADAG Equivalent |
|---------------------|------------------|
| `MAT::Node*` pointer | `NodeId` (size_t) |
| `MAT::Mutation` | `CompactGenome` entry + parent `EdgeMutations` entry |
| `Mutations_Collection` | `ContiguousMap<MutationPosition, ...>` |
| `node->mutations` | `edge.GetEdgeMutations()` on parent edge |
| `node->parent` | `dag.Get(node).GetParents()` (single parent in tree) |
| `node->children` | `dag.Get(node).GetChildren()` |
| `dfs_index` | Must precompute (not stored in MADAG) |
| `level` | Must precompute (not stored in MADAG) |
| `boundary1_all_major_allele` | New per-node data (AlleleCount) |
| Tree mutation | Overlay or mutable storage |
| `potential_crosses` | Vector indexed by NodeId or BFS index |

## File Organization Suggestion

```
include/larch/spr/
    move_enumerator.hpp        # Bounded search move enumeration
    allele_counts.hpp          # Per-node allele frequency tracking
    range_tree.hpp             # Range tree for pruning
    conflict_resolver.hpp      # Move conflict resolution
    working_tree.hpp           # Mutable working tree for compounding
    native_optimize.hpp        # Main optimization loop
include/larch/impl/spr/
    move_enumerator_impl.hpp   # Implementation
    allele_counts_impl.hpp
    range_tree_impl.hpp
    ...
test/
    test_move_enumerator.cpp   # Tests comparing with full Fitch
    test_native_optimize.cpp   # End-to-end optimization tests
```

## Testing Strategy

1. **Unit tests**: Each component tested against ground truth
   - Allele counts: verify against full tree scan
   - Range tree: verify pruning doesn't miss profitable moves
   - Incremental scoring: verify against `ParsimonyOnlyScoringBackend`

2. **Integration tests**: Compare optimization quality
   - Run native optimizer and matOptimize-based optimizer on same input
   - Native should find at least as good parsimony scores
   - Use the existing `sample_dag.hpp` test data

3. **Regression tests**: Ensure DAG integrity
   - All trees extractable from DAG are valid
   - Edge mutations consistent with compact genomes
   - Parsimony scores correct

4. **Performance benchmarks**: Measure speedup
   - Time per iteration vs matOptimize-based pipeline
   - Profile to identify bottlenecks
   - Compare on various tree sizes (100, 1k, 10k, 100k nodes)

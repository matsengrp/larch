# Existing Larch Infrastructure

What already exists in larch that can be reused for native MADAG optimization.

## SPR Move Evaluation

### ParsimonyOnlyScoringBackend

**Location**: `include/larch/spr/scoring_backend.hpp` (line 283+),
`include/larch/impl/spr/scoring_backend_impl.hpp` (line 330+)

Computes parsimony score change for an SPR move using native Fitch algorithm:

```cpp
// Usage:
auto spr = AddSPRStorageWithBackend<ParsimonyOnlyScoringBackend>(sampled_tree);
spr.View().InitHypotheticalTree(src_id, dst_id, lca_id);
int score_change = spr.View().GetScoreChange();
```

- Collects variable sites from edge mutations
- Runs Fitch on original and hypothetical trees
- Computes delta score
- Uses `uint8_t` one-hot Fitch sets per site per node
- **Limitation**: Full recomputation for each move (not incremental)

### SPR Overlay System

**Location**: `include/larch/spr/spr_view.hpp`, `include/larch/impl/spr/spr_view_impl.hpp`

Creates a non-destructive overlay representing the tree after an SPR move:

- **HypotheticalNode**: Tracks which nodes have changed topology or compact
  genomes. Methods: `IsMoveSource()`, `IsMoveTarget()`, `IsMoveNew()`,
  `IsLCAAncestor()`.
- **HypotheticalTree**: Wraps the scoring backend. `InitHypotheticalTree()`
  sets up the overlay and runs scoring.
- **Overlay storage**: Changes stored separately from the original tree. The
  original is never modified.

### Fragment Construction

**MakeFragment()** (spr_view_impl.hpp line 324+):
1. Identifies the oldest changed node's parent
2. Calls `PreorderComputeCompactGenome()` on changed nodes
3. Calls `CollapseEmptyFragmentEdges()` to remove non-mutating edges
4. Returns `FragmentStorage` with only the affected nodes and edges

**CollapseEmptyFragmentEdges()** (line 417+):
- Removes edges where parent and child have identical compact genomes
- Recursively finds the "final parent" of collapsed chains
- Deduplicates node and edge lists

## Random Move Generation

### RandomMoveGenerator

**Location**: `include/larch/spr/random_moves.hpp`

Generates valid (src, dst, lca) triples on a sampled tree:

```cpp
RandomMoveGenerator gen(sampled_tree.Const(), seed);
auto [src, dst, lca] = gen.GenerateMove();
```

- Validates moves: dst not in src's subtree, not src's parent, not root
- Supports distance constraints
- Uses `IsDescendant()` for validation (walks parent pointers)
- Retries up to `max_attempts` on invalid combinations

## Tree Sampling

### SubtreeWeight

**Location**: `include/larch/subtree/subtree_weight.hpp`

Samples trees from the DAG with optimal parsimony:

```cpp
SubtreeWeight<ParsimonyScore, MergeDAG> weight(dag);
auto tree = weight.MinWeightSampleTree({});
```

- Bottom-up DP computing optimal subtree weights
- Caches weights indexed by NodeId
- `MinWeightSampleTree()`: Samples a tree achieving minimum parsimony
- `MinWeightCount()`: Counts how many optimal trees exist

### ParsimonyScore

**Location**: `include/larch/subtree/parsimony_score.hpp`

Weight operations for parsimony:
- `ComputeEdge()`: Returns edge mutation count
- `WithinCladeAccumOptimum()`: Finds min-weight edges in a clade
- `BetweenClades()`: Sums weights across children
- `AboveNode()`: Adds edge weight to subtree weight

## Merge System

### Merge Class

**Location**: `include/larch/merge/merge.hpp`

Grows the DAG by incorporating new trees or fragments:

```cpp
Merge merge(reference_sequence);
merge.AddDAG(fragment.View());   // Single fragment
merge.AddDAGs(fragments);         // Batch of fragments
merge.ComputeResultEdgeMutations();
```

Key internals:
- `result_nodes_`: HashMap<NodeLabel, NodeId> -- deduplicates nodes
- `result_edges_`: HashMap<EdgeLabel, EdgeId> -- deduplicates edges
- **NodeLabel**: (CompactGenome, LeafSet, SampleId) -- unique node identity
- **EdgeLabel**: (parent NodeLabel, child NodeLabel, EdgeMutations)

Nodes with identical NodeLabels from different trees merge into one DAG node.
This is the core mechanism for DAG growth.

### Batched Merge

**Location**: Recent additions (commit 2b2c913)

`AddDAGs()` accepts multiple fragments at once for batch processing. This
avoids redundant recomputation between individual `AddDAG()` calls.

## Data Structures

### CompactGenome

**Location**: `include/larch/compact_genome.hpp`

Mutations from reference sequence:
```cpp
ContiguousMap<MutationPosition, MutationBase> data;
```

Key operations:
- `GetBase(pos, refseq)`: Base at position (from mutations or reference)
- `DifferingSites(other)`: Positions where two genomes differ
- `ToEdgeMutations(parent_cg)`: Derive edge mutations from parent-child CG pair
- `ApplyChanges(changes)`: Update genome with a set of changes

### EdgeMutations

**Location**: `include/larch/edge_mutations.hpp`

Parent-to-child mutations on an edge:
```cpp
ContiguousMap<MutationPosition, std::pair<MutationBase, MutationBase>> data;
```

Each entry: (position → (parent_base, child_base))

### ContiguousMap / ContiguousSet

**Location**: `include/larch/contiguous_map.hpp`, `contiguous_set.hpp`

Sorted vectors for cache-efficient lookup:
- O(log n) lookup via `std::lower_bound`
- O(n) merge/union operations
- Excellent cache locality for small-to-medium datasets
- Used throughout for mutations, genomes, and site tracking

### NodeId / EdgeId

**Location**: Strongly-typed `size_t` wrappers

Used to address nodes and edges in MADAG. More cache-friendly than pointers
(contiguous allocation) and naturally compatible with vector indexing.

## Optimization Loop Infrastructure

### OptimizeDAGWithRandomMoves

**Location**: `include/larch/spr/random_optimize.hpp`

Simple sequential optimization loop:
```
for each iteration:
    ComputeResultEdgeMutations()
    Sample min-weight tree
    Generate random moves
    For each move: create SPR overlay → make fragment → merge
    Merge sampled tree
```

### OptimizeDAGParallelRadius

**Location**: Same file

Three-phase parallel optimization with radius:
- Phase 1 (PARALLEL): Generate moves + SPR overlays
- Phase 2 (PARALLEL): Create fragments
- Phase 3 (BATCH): Merge all fragments

Uses exponential radius increases (1, 2, 4, 8, ...) similar to matOptimize.

## What's Missing for Efficient Enumeration

1. **Per-node allele count tracking**: Major allele + boundary1 per position
2. **Range tree / pruning structure**: For bounded search
3. **Incremental score update**: Without full Fitch recomputation
4. **DFS index precomputation**: For O(1) ancestor checks
5. **Conflict resolution**: For parallel move application
6. **Working tree mutation**: For compounding improvements within a radius
7. **Condensing**: Grouping identical leaf siblings (matOptimize does this)

## What Can Be Directly Reused

1. **Random move validation logic**: `IsDescendant()`, topology checks
2. **Fragment construction**: `MakeFragment()`, `CollapseEmptyFragmentEdges()`
3. **Merge system**: `Merge::AddDAG()`, `AddDAGs()`
4. **Tree sampling**: `SubtreeWeight::MinWeightSampleTree()`
5. **Data structures**: `ContiguousMap`, `CompactGenome`, `EdgeMutations`
6. **Parallel framework**: `ParallelForEach`, Taskflow integration
7. **Score verification**: `ParsimonyOnlyScoringBackend` for testing correctness

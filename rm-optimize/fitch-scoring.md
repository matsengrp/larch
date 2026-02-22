# Fitch Scoring and Parsimony

## Parsimony Scoring Basics

Parsimony score = total number of mutations across all edges in the tree.
Each edge where child has a different base than parent at some position
contributes +1 to the score.

The goal of optimization is to minimize this score.

## Fitch Algorithm

Bottom-up dynamic programming that finds optimal nucleotide assignments:

1. **Leaf nodes**: Fitch set = the observed nucleotide (e.g., {A})
2. **Internal nodes**:
   - Compute intersection of all children's Fitch sets
   - If intersection is non-empty: use it (no cost)
   - If intersection is empty: use union (cost +1 per pairwise disagreement)
3. **Top-down assignment**: Root picks any allele in its Fitch set. Each child
   picks from its Fitch set, preferring the parent's choice if possible.

## matOptimize's Incremental Fitch Tracking

matOptimize doesn't recompute Fitch from scratch for each candidate move.
Instead, it tracks how allele **counts** change incrementally.

### Key Concepts

#### Major Allele
The most frequent allele among a node's children at a given position.
Stored in `boundary1_all_major_allele` (lower 4 bits).

If the node adopts the major allele, it minimizes the number of children that
disagree (and thus minimizes parsimony cost at this node).

#### Boundary1 Allele
The second-most-frequent allele -- one that has count = (major_count - 1).
Stored in `boundary1_all_major_allele` (upper 4 bits).

Why track this? When a child is removed or added, the boundary1 allele might
become the new major allele (or vice versa), changing the parsimony cost.

#### Sensitive Increment
An allele whose count increase could improve parsimony. This happens when the
allele is already major or boundary1, meaning one more count could change which
allele the node adopts.

### State at Each Node (per position)

```
MAT::Mutation:
  position           - genomic site
  par_nuc             - parent's nucleotide (one-hot)
  mut_nuc             - this node's nucleotide (one-hot)
  all_major_allele    - major allele set among children (one-hot, may have ties)
  boundary1_allele    - alleles at count = major_count - 1 (one-hot)
```

### Incremental Update Rules

When **removing** a child with allele X from a node's children:

```
if X is the sole major allele:
    major_count decreases by 1
    if boundary1 is non-empty:
        boundary1 alleles become new major (tied with X)
    parsimony_change: may be +1 or 0
if X is one of tied major alleles:
    X drops to boundary1 or below
    remaining tied majors stay
    parsimony_change: 0 (still have other major choices)
if X is boundary1:
    X drops below boundary1
    parsimony_change: 0
if X is below boundary1:
    no effect on scoring
```

When **adding** a child with allele X:

```
if X is already major:
    major_count increases, no parsimony change
if X is boundary1:
    X ties for major, may reduce parsimony if X matches parent
if X is below boundary1:
    may become boundary1, no immediate parsimony change
if X is not present:
    count starts at 1, usually below boundary1
```

## Larch's Native Fitch Implementation

Larch has its own Fitch implementation in `ParsimonyOnlyScoringBackend`
(`include/larch/spr/scoring_backend.hpp` and `impl/spr/scoring_backend_impl.hpp`).

### How It Works

1. **Collect variable sites**: All positions with mutations in the tree
2. **Pre-move Fitch pass**: Bottom-up on original tree, storing Fitch sets
3. **Post-move Fitch pass**: Bottom-up on hypothetical (overlay) tree
4. **Score comparison**: `score_change = new_score - old_score`

### Fitch Set Storage

Per-node, per-site `uint8_t` with one-hot encoding:
```
A = 0x1, C = 0x2, G = 0x4, T = 0x8
```

Bitwise operations for Fitch:
```cpp
intersection = child1_set & child2_set;
union_set = child1_set | child2_set;
if (intersection != 0)
    node_set = intersection;  // No cost
else
    node_set = union_set;     // Cost +1
```

### Limitations for Optimization

Larch's current Fitch implementation:

- **Recomputes from scratch** for each candidate move. Fine for random moves
  (only evaluating a few), but too expensive for exhaustive enumeration.
- **No boundary1/major allele tracking**: Cannot do incremental updates.
- **No range tree support**: Cannot prune unpromising subtrees.

### What's Needed for Efficient Enumeration

To match matOptimize's performance, larch needs:

1. **Per-node allele counts** at each variable site (not just Fitch sets)
2. **Major allele + boundary1 tracking** for incremental updates
3. **Incremental score change computation** without full tree traversal
4. **Range tree construction** from MADAG's edge mutations

## Scoring a Candidate Move: The Three Components

For an SPR move (src, dst, LCA):

### Component 1: Removal Score

Removing src from its parent changes the parent's allele counts. This propagates
upward: the parent's major allele may change, affecting the grandparent, etc.

Computed during the upward traversal in move enumeration.

### Component 2: Split Score

Inserting src as sibling of dst creates a new internal node. The new node has
two children (src and dst). Its Fitch set = intersection of src's and dst's
alleles (or union if they disagree).

The cost depends on:
- Whether src's allele matches dst's allele at each position
- Whether the new node's chosen allele matches dst's parent's allele
- Whether this changes dst's parent's major allele set

### Component 3: Above-LCA Score

Changes at the LCA propagate upward. The LCA loses one child branch (src's
old position) and gains a modified subtree. Its major allele may change,
propagating the effect toward the root.

**Total score change = removal + split + above_LCA**

Negative total = improvement. Zero = neutral (accepted in drift mode).

## Comparison: matOptimize vs Larch Fitch

| Aspect | matOptimize | Larch Native |
|--------|-------------|--------------|
| State tracking | Major allele + boundary1 per node | Fitch sets per node |
| Update strategy | Incremental (O(changed sites)) | Full recomputation (O(all variable sites)) |
| Pruning | Range tree + lower bounds | None |
| Data structure | Mutations_Collection (sorted vec) | ContiguousMap + uint8_t vectors |
| Parallelism | TBB flow graph | ParallelForEach / Taskflow |
| Suitable for | Exhaustive enumeration | Few random moves |

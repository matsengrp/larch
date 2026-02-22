# Move Enumeration: The Bounded Search Algorithm

This is the most complex and performance-critical part of matOptimize. It
answers: "Given a source node, which destination nodes yield profitable SPR
moves?"

## Naive Approach (What Random Moves Do)

Pick random (src, dst) pairs, compute full parsimony for each. This is O(N) per
move evaluation and requires trying O(N) destinations. Total: O(N^2) per source
node -- prohibitive for large trees.

## matOptimize's Approach: Incremental Scoring with Bounded Pruning

Instead of computing full parsimony for each candidate, matOptimize:

1. **Incrementally** tracks how allele counts change as it traverses the tree
2. **Prunes** subtrees that provably cannot contain profitable destinations
3. Achieves sub-quadratic search in practice (though worst-case is still O(N^2))

## Core Data Structure: Mutation_Count_Change

Tracks how removing/adding a subtree changes allele frequencies at each genomic
position:

```cpp
struct Mutation_Count_Change {
    int position;           // Genomic site
    nuc_one_hot decremented_allele;  // Allele losing a count (4-bit one-hot)
    nuc_one_hot incremented_allele;  // Allele gaining a count (4-bit one-hot)
    nuc_one_hot par_state;           // Parent's current allele
};
```

One-hot encoding: A=0x1, C=0x2, G=0x4, T=0x8. Bitwise ops enable fast set
operations.

## The Two-Phase Traversal

### Phase 1: Upward (src toward root)

Starting at the source node, walk toward the root. At each ancestor:

1. **Remove src's contribution**: The ancestor loses one child's allele at each
   mutated position. Update `Mutation_Count_Change` accordingly.
2. **Check this node as LCA**: If we prune src and regraft it somewhere in the
   sibling subtrees of this ancestor, this ancestor becomes the LCA. Evaluate
   profitability.
3. **Search sibling subtrees** (Phase 2) for profitable destinations.
4. **Propagate upward**: Merge mutation changes to prepare for the next ancestor
   level.

At each step, the algorithm maintains a running set of `Mutation_Count_Change`
objects representing "what happens to parsimony at this level if src is removed."

### Phase 2: Downward (from LCA into sibling subtrees)

For each potential LCA found in Phase 1, search its non-src children downward:

1. At each node in the subtree, compute the **split cost**: the parsimony change
   from inserting a new internal node that makes src a sibling of this node.
2. Use the **range tree** to determine if going deeper can possibly improve the
   score.
3. If the lower bound already exceeds the best known score, prune.

## Range Tree: The Key Pruning Structure

### What it stores

For each genomic position, a hierarchical index over the tree's DFS ordering:

```cpp
struct range_tree_node {
    uint32_t dfs_start_idx, dfs_end_idx;  // DFS range covered
    LEVEL_T min_level[4];                  // Min depth per allele (A,C,G,T)
    IDX_TREE_IDX_T children_start_idx;
    IDX_TREE_IDX_T parent_idx;
};
```

### How pruning works

For a mutation at position P where we're incrementing allele X:
- Look up `min_level[X]` in the range tree node covering the current subtree
- If `min_level[X] > radius_left`, no node within remaining radius has allele X
  at position P. This mutation **cannot** improve parsimony in this subtree.
- Increment the lower bound counter. If lower_bound > best_score, prune entirely.

### Why it works

The range tree encodes: "What is the shallowest tree depth at which I can find a
node that has allele X at position P?" If all such nodes are deeper than the
remaining search radius allows, there is no point searching this subtree for
position P's contribution.

### Construction

Built once per search iteration from the tree's DFS-ordered nodes:
- For each genomic position, collect all nodes with a mutation
- Build hierarchical intervals over DFS ranges
- Precompute `min_level` per allele at each interval

This is O(N * S) where S = number of variable sites (positions with mutations).

## Lower Bound Tracking

The extended structures carry lower bounds:

```cpp
struct Mutation_Count_Change_W_Lower_Bound_Downward {
    // ... base fields ...
    nuc_one_hot par_sensitive_increment;  // Alleles that CAN reduce parsimony
    LEVEL_T next_level;                    // Deepest useful level
    IDX_TREE_IDX_T idx;                    // Position in range tree
};
```

**Sensitive increment**: An allele is "sensitive" if incrementing it could reduce
parsimony. This happens when the allele is already in the major or boundary1
set. If incrementing a non-sensitive allele, it cannot help.

The lower bound accumulates across all positions. If the cumulative lower bound
exceeds the best known improvement, the entire subtree is pruned.

## Score Computation at a Candidate Destination

When a candidate dst is reached, the final score combines three components:

### 1. Score change from removing src from its current location

Computed during the upward phase. Removing src changes its parent's major allele
set, which may propagate upward.

### 2. Score change from the split at dst

Inserting src as sibling of dst creates a new internal node. The cost depends on
whether src's and dst's alleles agree at each position:

```
if src_allele == dst_allele:
    split_cost = 0  (new node inherits shared allele)
elif src_allele == parent_allele:
    split_cost = 0  (new node uses parent allele, dst already has mutation)
else:
    split_cost = +1  (new node must choose, one child pays a mutation)
```

### 3. Score change above LCA

Changes to the LCA's major allele set propagate upward. The algorithm walks from
LCA to root, checking at each ancestor if the major allele changed.

The functions `check_move_profitable_dst_not_LCA()` and
`check_move_profitable_LCA()` combine these three components.

## DFS Ordering for Fast Ancestor Checks

Every node has precomputed:
- `dfs_index`: Position in preorder traversal
- `dfs_end_index`: End of subtree in preorder

This enables O(1) ancestor test:
```
is_ancestor(a, b) = (b.dfs_index >= a.dfs_index) && (b.dfs_index <= a.dfs_end_index)
```

Used heavily during move enumeration to avoid walking parent pointers.

## Handling Internal vs Leaf Nodes

The parsimony change differs:

**Leaf nodes** (terminal):
```
if incremented_allele not in par_state: +1  (new disagreement with parent)
if decremented_allele was only in par_state: -1  (removed disagreement)
```

**Internal nodes**:
```
if incremented_allele in par_state: -1  (child can now follow parent)
if decremented_allele in par_state: +1  (child can no longer follow parent)
```

The distinction: leaves count the mutation itself (does child differ from
parent?), while internal nodes count whether children can agree with the
parent's choice.

## Complexity Analysis

- **Worst case per source**: O(N * S) where N = nodes, S = variable sites
- **With pruning**: Much better in practice. The range tree prunes large subtrees
  early, and the lower bound prevents exploring dead ends.
- **Total per iteration**: O(K * average_search_cost) where K = nodes to search
- **Range tree construction**: O(N * S) one-time per iteration

The bounded search is what makes matOptimize practical for trees with 100k+
nodes. Without it, the algorithm would be quadratic and infeasible.

## Mapping to MADAG

For a native MADAG implementation, the key questions are:

1. **Mutation representation**: MADAG uses `EdgeMutations` (parent_base,
   child_base pairs) and `CompactGenome` (mutations from reference). The
   `Mutation_Count_Change` tracking needs to work with these.

2. **DFS ordering**: MADAG doesn't precompute DFS indices. Either precompute
   them on each sampled tree, or find alternatives.

3. **Range tree**: Can be built from MADAG's edge mutations. Each edge's
   mutations give the positions and alleles needed.

4. **Fitch sets**: Larch already has native Fitch computation
   (`ParsimonyOnlyScoringBackend`). The incremental version (tracking boundary1
   and major alleles) would need to be added.

5. **Node identity**: MADAG nodes have NodeId (size_t). matOptimize uses
   pointers. NodeId is simpler and cache-friendlier.

See [implementation-plan.md](implementation-plan.md) for the mapping in detail.

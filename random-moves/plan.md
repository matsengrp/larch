# Pure-MADAG Random Move Optimization Loop

## Context

The current optimization loop (`optimize_dag_direct` in `produce_mat_impl.hpp`)
requires converting between larch's MADAG and matOptimize's MAT::Tree at every
step: BuildMAT before optimization, Move_Found_Callback with MAT::Node*
pointers during, BuildFromMAT after. The MATConversion feature, MATView layer,
condensation cycle, and Fitch set plumbing through the scoring backend all exist
to bridge this gap.

The goal is to implement the same sample→move→fragment→merge loop using only
larch-native MADAG structures and NodeIds, with a random move engine (like the
mock optimizer) that doesn't depend on matOptimize at all. This also requires
refactoring the HypotheticalTree/scoring backend layer to cleanly support
MAT-free operation.

## Approach

### Part 1: New scoring backend — `ParsimonyOnlyScoringBackend`

**File: `include/larch/spr/scoring_backend.hpp`** (add new class)
**File: `include/larch/impl/spr/scoring_backend_impl.hpp`** (add impl)

A minimal backend that uses larch's own `SankoffScorer` (from
`include/larch/subtree/sankoff.hpp`) to compute Fitch sets and assign optimal
bases at internal nodes after an SPR move. No MAT dependency.

Interface (mirrors `MatOptimizeScoringBackend` shape):
- `Initialize(dag, src, dst, lca)` — stores NodeIds, runs Sankoff on the
  affected subtree to compute per-node optimal base sets
- `GetMoveSource/Target/LCA()` — return NodeIds
- `GetScoreChange()` — return parsimony delta computed by Sankoff
- `GetFitchSetAtSite(dag, node, site, ...)` — return `FitchSet` from Sankoff DP
  table (the set of equally-optimal bases at that site for that node)
- `GetFitchSetParts(dag, node, ...)` — return `{empty Mutations_Collection,
  optional changes map}`. The changes map can be synthesized from the Sankoff DP
  costs (sites where the optimal set changed from the pre-move state).
- `GetSitesWithScoringChanges(node)` — sites where Sankoff changed the optimal set
- `HasScoringChanges(node)` — whether any sites changed
- `SelectBase(fitch_set, old_base, parent_base)` — same logic as matOptimize
  backend (prefer parent, then old, then first available)

Key design: Run `SankoffScorer` on the hypothetical tree (post-ApplyMove overlay)
to get optimal bases. The scorer works on any DAG view — the overlay created by
`ApplyMove` is a valid tree view.

The Sankoff pass is a full bottom-up + top-down over the tree. For efficiency,
we could restrict it to the affected subtree (nodes between the fragment root
and leaves). But for random moves the full-tree pass is acceptable as a starting
point — optimize later if needed.

### Part 2: Refactor `spr_view_impl.hpp` for MAT-free backends

**File: `include/larch/impl/spr/spr_view_impl.hpp`**

Three changes to remove hard MAT dependencies:

1. **`IsMATRoot()` (line 2-8)**: Currently calls `GetMATNode()`. Replace with a
   topology check that works for all backends:
   ```cpp
   bool IsMATRoot() const {
     auto& node = static_cast<const CRTP&>(*this);
     auto old = node.GetOld();
     if (old.IsUA()) return false;
     return old.GetSingleParent().GetParent().IsUA();
   }
   ```
   This checks: "is this node the tree root?" (i.e., its parent is the UA).
   Rename to `IsTreeRoot()` for clarity; keep `IsMATRoot` as alias.

2. **`ComputeNewCompactGenome()` (line 125)**: Remove `Assert(node.HaveMATNode())`.
   This assert is wrong for MAT-free backends. The function body works fine
   without MAT — it delegates to backend's `GetFitchSetAtSite` and `SelectBase`,
   which our new backend provides.

3. **`HypotheticalTree::Data` NodeId constructor (line 1233-1252)**: Currently has
   `if constexpr (is_same<Backend, MLScoringBackend>)` to decide whether to call
   MAT-based or ML-based init. Add `ParsimonyOnlyScoringBackend` to the non-MAT
   branch:
   ```cpp
   if constexpr (std::is_same_v<Backend, MLScoringBackend<DAG>> ||
                 std::is_same_v<Backend, ParsimonyOnlyScoringBackend<DAG>>) {
     backend_.Initialize(dag, src, dst, lca);
   } else {
     // existing MatOptimize path with GetMATNode conversion
   }
   ```

**File: `include/larch/spr/spr_view.hpp`**

4. **`GetChangedFitchSetMap()` (line 151-155)**: Returns
   `ContiguousMap<MATNodePtr, ...>`. For MAT-free backends this should return an
   empty map. The existing `MatOptimizeScoringBackend` has this method. Add a
   stub to `ParsimonyOnlyScoringBackend` that returns a static empty map.

5. **`Profitable_Moves move_` field (line 101)**: Stays as-is — it's
   default-initialized (nullptr members) for non-matOptimize paths, and nothing
   in the MAT-free path reads it.

### Part 3: Random move generator

**New file: `include/larch/spr/random_moves.hpp`**

```cpp
template <typename DAG>
class RandomMoveGenerator {
public:
  struct Move { NodeId src, dst, lca; };

  RandomMoveGenerator(DAG dag, std::optional<uint32_t> seed = std::nullopt);
  std::optional<Move> GenerateMove();

private:
  DAG dag_;
  std::vector<NodeId> searchable_nodes_;  // non-UA, non-root
  std::mt19937 gen_;
};
```

Logic adapted from mock `optimize_inner_loop` in `optimize.hpp:458-554`, but
using NodeIds and larch DAG traversal:

1. Collect all nodes that are not UA into `searchable_nodes_`
2. Pick random `src` (must have a parent, must not be root)
3. Pick random `dst` (must not be src, not parent of src, not descendant of src)
4. Descendant check: walk `dst`'s ancestors via `GetSingleParent().GetParent()`
   until UA, checking none equal `src`
5. Find LCA via `FindLCA(dag.Get(src), dag.Get(dst))` from `lca.hpp`
6. Validate: LCA parent must not be UA (the mock optimizer has this check too)
7. Return `{src, dst, lca.lca}`

### Part 4: The pure-MADAG optimization loop

**New file: `include/larch/spr/random_optimize.hpp`**

```cpp
template <typename MergeT>
void OptimizeDAGWithRandomMoves(
    MergeT& merge,
    size_t num_iterations,
    size_t moves_per_iteration,
    bool collapse_empty_fragment_edges = true,
    std::optional<uint32_t> seed = std::nullopt);
```

Each iteration:
1. `merge.ComputeResultEdgeMutations()`
2. Sample a tree: `SubtreeWeight<ParsimonyScore, ...> weight{merge.GetResult()}`
   → `auto sampled = weight.MinWeightSampleTree({})` (returns
   `SampledDAGStorage<MADAGStorage<>>` with `MappedNodes`, no MATConversion)
3. `sampled.View().RecomputeCompactGenomes(true)`
4. `sampled.View().SampleIdsFromCG()`
5. Create `RandomMoveGenerator` on `sampled.View()`
6. For each move:
   a. `auto spr = AddSPRStorageWithBackend<ParsimonyOnlyScoringBackend<...>>(sampled.View())`
   b. `spr.View().InitHypotheticalTree(move.src, move.dst, move.lca)`
   c. If init succeeded: `auto fragment = spr.View().MakeFragment()`
   d. Collect fragment
7. `merge.AddDAGs(fragments)` — batch-merge all fragments
8. Also merge the sampled tree: `merge.AddDAG(sampled.View())`

No MAT::Tree, no MATConversion, no BuildMAT/BuildFromMAT, no
Move_Found_Callback, no BatchingCallback.

### Part 5: Test

**New file: `test/test_random_optimize.cpp`**

Using the sample DAG from `test/sample_dag.hpp`:
- **Test 1**: Create DAG, run `OptimizeDAGWithRandomMoves` for a few iterations,
  verify DAG size grows (more trees merged in)
- **Test 2**: Verify merged DAG is valid — all edges have valid mutations,
  CompactGenomes are consistent
- **Test 3**: Verify parsimony score is computable on result trees
- **Test 4**: Single random move — generate one, apply via SPR, verify fragment
  is valid

Register tests with tags `{"random-moves"}`.

## Files Summary

| File | Action | Description |
|------|--------|-------------|
| `include/larch/spr/scoring_backend.hpp` | Modify | Add `ParsimonyOnlyScoringBackend` class |
| `include/larch/impl/spr/scoring_backend_impl.hpp` | Modify | Add impl for new backend |
| `include/larch/impl/spr/spr_view_impl.hpp` | Modify | Fix `IsMATRoot`, remove `HaveMATNode` assert, update Data constructor `if constexpr` |
| `include/larch/spr/random_moves.hpp` | Create | Random move generator |
| `include/larch/spr/random_optimize.hpp` | Create | Pure-MADAG optimization loop |
| `test/test_random_optimize.cpp` | Create | Tests |

## Existing code reused

- `SankoffScorer` from `include/larch/subtree/sankoff.hpp` — Fitch/Sankoff for
  base assignment
- `FindLCA` from `include/larch/spr/lca.hpp` — LCA computation
- `ApplyMove` from `spr_view_impl.hpp` — already NodeId-based
- `MakeFragment` from `spr_view_impl.hpp` — already NodeId-based
- `SubtreeWeight::MinWeightSampleTree` from
  `include/larch/subtree/subtree_weight.hpp` — tree sampling
- `AddSPRStorageWithBackend<>` from `spr_view.hpp:241` — parameterized SPR storage
- `InitHypotheticalTree(NodeId, NodeId, NodeId)` from `spr_view_impl.hpp:1158`
  — the ML/NodeId init path

## Verification

Build and test in the `build-deb-mock-asan` directory (USE_USHER=no,
USE_ASAN=yes):

```bash
cd build-deb-mock-asan
make -j32 larch-test
./bin/larch-test "random-moves:.*"
```

This configuration uses mock MAT types from `optimize.hpp` (no real matOptimize
compiled), which confirms the loop works without matOptimize. ASan catches any
memory issues.

Also verify existing tests still pass:
```bash
./bin/larch-test "SPR:.*"
./bin/larch-test "Sankoff:.*"
```

# Implementation Notes

Quirks, gotchas, and discoveries found during implementation.

## Overlay System: Mutable vs Const Views

The overlay system's `is_mutable` flag determines write behavior for non-overlaid nodes:

```cpp
constexpr bool is_mutable = (not std::is_const_v<OverlayStorageType>) and
                            OverlayStorageType::TargetView::is_mutable;
```

- **Mutable overlay + non-overlaid node write** → `Fail("Can't modify non-overlaid node")`
- **Const overlay + non-overlaid node read** → falls through to target storage

**Key pattern**: Use `dag.Const()` when reading through an overlay to avoid triggering the mutable write check. This is used in `ParsimonyOnlyScoringBackend::Initialize()` when running Fitch on the overlay tree.

## MappedNodes and Fragment Merging

`SampledDAGStorage` (from `MinWeightSampleTree`) extends with `MappedNodes`. This propagates through the SPR overlay chain to fragments:

```
SampledDAGStorage (has MappedNodes)
  → SPRStorage overlay
    → FragmentStorage (delegates to SPR overlay)
```

When `merge.AddDAG(fragment.View())` is called, `MergeNodes` checks `contains_element_feature<Node, MappedNodes>` at compile time and calls `SetOriginalId()`. For non-overlaid fragment nodes, this triggers the "Can't modify non-overlaid node" error.

**Fix**: Before merging, overlay `MappedNodes` on all non-appended SPR nodes that the fragment references:

```cpp
auto fragment = spr.View().MakeFragment();
for (auto node : fragment.View().GetNodes()) {
  auto spr_node = spr.View().Get(node.GetId());
  if (not spr_node.IsAppended()) {
    spr_node.template SetOverlay<MappedNodes>();
  }
}
merge.AddDAG(fragment.View());
```

The batching callback avoids this because its SPR chain uses `MATConversionStorage<MergeDAGStorage<>>` which does NOT include `MappedNodes`.

## Backend Dispatch via C++20 `requires`

The SPR `Data` constructor uses `if constexpr` to select between matOptimize-based and pure-MADAG backends. Using `is_same_v` fails because the backend template parameter's DAG type is a storage type while the actual DAG is a view type.

**Fix**: Use a C++20 `requires` expression to detect the backend by its `Initialize` signature:

```cpp
if constexpr (requires(Backend b, DAGView d, Profitable_Moves m,
                       std::vector<Node_With_Major_Allele_Set_Change> c) {
                b.Initialize(d, m, c);
              }) {
  // matOptimize backend path
} else {
  // ParsimonyOnly backend path
}
```

## IsMATRoot Without MAT Dependency

The original `IsMATRoot()` called `GetMATNode()` which requires matOptimize. Replaced with a topology-based check:

```cpp
bool IsMATRoot() const {
  auto old = node.GetOld();
  if (old.IsUA()) { return false; }
  return old.GetSingleParent().GetParent().IsUA();
}
```

## Leaf CGs in Merge Results

Leaf node CompactGenomes are intentionally **empty** in the merge result DAG. Leaf sequences are identified solely by SampleId, with authoritative CGs stored in `sample_id_to_cg_map`. The `ComputeResultEdgeMutations` function uses this map for leaf children, not the (empty) node CGs.

This means edge mutations for leaf children are consistent with the authoritative CGs, NOT with what `node.GetCompactGenome().GetBase()` returns (which falls back to the reference sequence for empty CGs).

## Fitch Algorithm Details

`ParsimonyOnlyScoringBackend` implements unweighted parsimony via Fitch:

1. **Variable sites**: Collected from original tree's edge mutations
2. **Leaf Fitch sets**: Singleton containing the leaf's base (one-hot encoded via `base_to_singleton`)
3. **Internal node sets**: Intersection of children if non-empty, else union (with parsimony score increment)
4. **Changed sites**: Compared between pre-move and post-move Fitch sets per node

The `GetFitchSetAtSite` method returns pre-computed post-move Fitch sets for internal nodes, and the original (pre-move) singleton for leaf nodes (leaves don't change in an SPR move).

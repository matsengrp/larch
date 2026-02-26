# matOptimize / Mock Optimizer Usage Outside Core Integration Files

References to matOptimize types and functions outside of `mat_conversion.hpp`,
`mat_view.hpp`, `dag_loader*`, and `usher_glue.hpp`.

---

## 1. The Optimization Loop: `produce_mat_impl.hpp`

**File:** `include/larch/impl/produce_mat_impl.hpp`

This is the implementation of `optimize_dag_direct()` — the main optimization
entry point declared in `usher_glue.hpp`. It orchestrates the full
sample-optimize-merge cycle.

**Usage:**
- `dag.GetMutableMAT()` — gets the `MAT::Tree&` from the MATConversion feature
- `save_mutation_annotated_tree(tree, ...)` — writes `.pb` snapshots before/after optimization
- `Original_State_t origin_states` — builds the initial state map
- `check_samples(tree.root, origin_states, &tree)` — validates leaf mutations
- `reassign_states(tree, origin_states)` — re-optimizes internal node labels
- `tree.condense_leaves(...)` / `tree.uncondense_leaves()` — condensation round-trip
- `optimize_inner_loop(all_nodes, tree, radius, ...)` — the actual matOptimize call (or mock random moves)
- `AddMATConversion(MADAGStorage<>::EmptyDefault())` — creates result storage
- `result.first.View().BuildFromMAT(result.second, ...)` — converts optimized tree back to DAG

---

## 2. SPR Scoring Backend: `scoring_backend.hpp` / `scoring_backend_impl.hpp`

**Files:** `include/larch/spr/scoring_backend.hpp`, `include/larch/impl/spr/scoring_backend_impl.hpp`

`MatOptimizeScoringBackend` translates matOptimize's move representation into
larch's NodeId-based system for hypothetical tree construction.

**Usage:**
- `Initialize(dag, Profitable_Moves&, vector<Node_With_Major_Allele_Set_Change>&)` — converts `MAT::Node*` pointers to `NodeId` via `GetNodeFromMAT()`
- `Mutation_Count_Change` — stored in `ContiguousMap<MutationPosition, Mutation_Count_Change>` per node, representing Fitch set changes
- `MATNodePtr`-keyed map (`mat_keyed_fitch_set_map_`) — legacy access path alongside the `NodeId`-keyed map
- `MAT::Mutations_Collection` — returned by `GetFitchSetParts()` for nodes with changed Fitch sets
- `node.GetMATNode()->mutations` — reads mutations directly from the MAT node

---

## 3. SPR Hypothetical Tree: `spr_view.hpp` / `spr_view_impl.hpp`

**Files:** `include/larch/spr/spr_view.hpp`, `include/larch/impl/spr/spr_view_impl.hpp`

`HypotheticalTree` represents what a tree would look like after applying an SPR
move, without actually modifying the DAG.

**Usage:**
- `HypotheticalTree::Data` stores a `Profitable_Moves move_` field (src/dst/LCA as `MAT::Node*`)
- Constructor takes `const Profitable_Moves&` and `vector<Node_With_Major_Allele_Set_Change>&`
- `InitHypotheticalTree(Profitable_Moves&, ...)` — initializes the hypothetical tree from the move
- `node.GetMATNode()->parent == nullptr` — checks if a node is root via its MAT pointer
- `move_.src = dag.Get(src).GetMATNode()` — resolves NodeId back to MAT pointer for the move
- `ContiguousMap<MATNodePtr, ContiguousMap<MutationPosition, Mutation_Count_Change>>` — Fitch set changes keyed by MAT node pointer

---

## 4. Batching Callback: `batching_callback.hpp` / `batching_callback_impl.hpp`

**Files:** `include/larch/spr/batching_callback.hpp`, `include/larch/impl/spr/batching_callback_impl.hpp`

`BatchingCallback<CRTP>` inherits from `Move_Found_Callback` and adds batching
of accepted moves into DAG fragments for merging.

**Usage:**
- Inherits `Move_Found_Callback` — matOptimize's callback interface
- `operator()(Profitable_Moves&, int, vector<Node_With_Major_Allele_Set_Change>&)` — accepts/rejects individual moves
- `operator()(MAT::Tree&)` — called at end of each radius with the full tree
- `OnReassignedStates(MAT::Tree&)` — called after `reassign_states`
- `CopyNode(const MAT::Node*, MAT::Node*)` / `CopyTree(const MAT::Tree&)` — deep-copies MAT trees for per-callback snapshots
- `MAT::Tree sample_mat_tree_` — stores a copy of the sampled MAT tree
- `SetSample(MAT::Tree&, string)` — copies a MAT tree into the callback
- `CreateMATViewStorage(MAT::Tree&, ...)` / `CreateMATStorage(MAT::Tree&, ...)` — builds MATConversion storage from a MAT tree
- `AddMATConversion(...)` — wraps merge/sample storage with MATConversion
- `node.GetMATNode()->node_id` — reads MAT node IDs for sample ID lookup
- `tree.condensed_nodes.at(node.GetMATNode()->node_id)` — resolves condensed node names

---

## 5. VCF Original State: `vcf_original_state.hpp`

**File:** `include/larch/vcf/vcf_original_state.hpp`

Builds the `Original_State_t` map from VCF input, used to initialize
`reassign_states`.

**Usage:**
- `Original_State_t MakeOriginalState(VcfFile&)` — return type is matOptimize's state map
- `MAT::get_nuc_id(char)` — converts IUPAC characters to matOptimize's nucleotide encoding
- `MAT::Mutation(chromosome, pos, ...)` — constructs matOptimize mutation objects
- `Mutation_Set` (from `Original_State_t`) — per-node set of mutations

---

## 6. `tools/larch-usher.cpp` — Main CLI Executable

The largest consumer of matOptimize types outside the library headers.

### Callback classes (all inherit `BatchingCallback`):

- **`Treebased_Move_Found_Callback`** (line 136): Accepts moves based on score
  thresholds, builds hypothetical trees, and merges accepted fragments. Uses
  `GetMATNode()->node_id` and `GetMAT().get_node(nid)` to cross-reference
  nodes between hypothetical trees and the mapped storage.

- **`Merge_All_Moves_Found_Callback`** (line 245): Merges all moves regardless
  of score (no hypothetical tree check).

- **`Merge_All_Profitable_Moves_Found_Callback`** (line 267): Merges only
  score-improving moves, with hypothetical tree verification. Stores
  `atomic<MAT::Tree*> sample_mat_` for concurrent access.

- **`Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback`** (line 389): Like
  above but keeps the sampled tree fixed across radii. Also stores
  `atomic<MAT::Tree*> sample_mat_`.

### Optimization loop (lines 930-994):

- `AddMATConversion(weight.SampleTree(...))` — wraps sampled trees with MATConversion
- `MAT::Tree mat` + `sample_tree.View().BuildMAT(mat)` — converts sampled DAG to MAT
- `optimize_dag_direct(sample_tree.View(), callback, callback, callback)` — runs the optimization loop with one of the callback classes above

---

## 7. Test Files

### `test/test_spr.cpp`
- `Empty_Callback : Move_Found_Callback` — no-op callback for basic SPR tests
- `Test_Move_Found_Callback : BatchingCallback<...>` — exercises the batching path
- `Single_Move_Callback_With_Hypothetical_Tree : Move_Found_Callback` — tests hypothetical tree construction from `Profitable_Moves`; uses `AddMATConversion`, `BuildFromMAT`, `GetNodeFromMAT`, `GetMAT`
- Calls `optimize_dag_direct` with both callback types

### `test/test_spr_after_move.cpp`
- `Test_Move_Found_Callback : Move_Found_Callback` — tests DAG state after SPR moves
- Uses `MATConversionStorage` typedef, `BuildFromMAT`, `reassign_states` callback
- Stores `atomic<MAT::Tree*> sample_mat_` for concurrent move handling
- `AddMATConversion(make_sample_dag())` + `BuildMAT(tree)` in test setup

### `test/test_matOptimize.cpp`
- `Test_Move_Found_Callback : Move_Found_Callback` — simple callback returning `map<MAT::Node*, CompactGenome>`
- `Larch_Move_Found_Callback : Move_Found_Callback` — full callback with `MATNodePtr` to `NodeId` mapping, `OnReassignedStates(MAT::Tree&)`, and hypothetical tree verification
- Tests `optimize_dag_direct` with radius callbacks that call `BuildFromMAT` and remap `MATNodePtr` to `NodeId`

### `test/test_hypothetical_tree_pseudocode.cpp`
- `Single_Move_Callback_With_Hypothetical_Tree : Move_Found_Callback` — tests hypothetical tree construction
- `AddMATConversion`, `BuildMAT`, `optimize_dag_direct`

### `test/test_mat_conversion.cpp`
- Tests `BuildMAT` and `BuildFromMAT` round-trip correctness
- Verifies sample IDs survive MAT conversion

### `test/test_mat_view.cpp`
- Tests MAT view layer: `BuildMAT`, `check_samples`, `reassign_states`, `condense_leaves`
- Verifies condensed and uncondensed representations

### `test/test_fileio_protobuf.cpp`
- Tests protobuf I/O by round-tripping through `AddMATConversion` + `BuildMAT` + `BuildFromMAT`

### `test/test_overlay.cpp`
- `AddMATConversion(dag)` + `BuildMAT(mat)` + `mat.condense_leaves(...)` — tests overlay behavior on condensed trees

### `test/test_vcf.cpp`
- `AddMATConversion(MakeSampleDAG())` + `BuildMAT(tree)` + `reassign_states(tree, orig_state)` — tests VCF-to-MAT pipeline

### `test/test_ambiguous_vcf.cpp`
- `Original_State_t` parameter in helper functions
- `MAT::Mutation` construction for comparing expected vs actual mutations
- `AddMATConversion` + `BuildMAT` in test setup

---

## Summary of matOptimize Types by Usage Category

| Type | Where Used | Purpose |
|------|-----------|---------|
| `MAT::Tree` | produce_mat, batching_callback, larch-usher, tests | The tree matOptimize operates on |
| `MAT::Node` / `MATNodePtr` | scoring_backend, spr_view, batching_callback, larch-usher | Cross-referencing nodes between MAT and DAG |
| `MAT::Mutation` / `Mutations_Collection` | scoring_backend, vcf_original_state, test_ambiguous_vcf | Per-node mutation data |
| `Profitable_Moves` | spr_view, scoring_backend, batching_callback, all callbacks | SPR move descriptor (src/dst/LCA) |
| `Move_Found_Callback` | batching_callback, all test/tool callbacks | Callback interface for move discovery |
| `Node_With_Major_Allele_Set_Change` | spr_view, scoring_backend, batching_callback, all callbacks | Per-node Fitch set changes from a move |
| `Mutation_Count_Change` | spr_view, scoring_backend | Per-site allele count delta |
| `Original_State_t` | produce_mat, vcf_original_state, test_mat_view, test_vcf | Initial state for reassign_states |
| `MATConversion` (larch Feature) | produce_mat, batching_callback, larch-usher, most tests | Storage feature linking DAG nodes to MAT nodes |
| `AddMATConversion()` | batching_callback, larch-usher, most tests | Wraps a DAG storage with MATConversion |

# Larch

Tools for inference and manipulation of history DAGs (phylogenetic trees).

## Code Style

- Use clang-format (config at `.clang-format`)
- Types and global/member functions: `CamelCase`
- Member data: `trailing_underscore_`

## Building

Multiple build directories for different compile-time option combinations. Naming convention:

- Start with `build-`
- `deb-` or `rel-` for Debug or RelWithDebInfo
- `mock-` if USE_USHER=no
- `asan-` or `tsan-` if USE_ASAN=yes or USE_TSAN=yes
- `serial-` if DISABLE_PARALLELISM=yes
- `asserts-` if KEEP_ASSERTS=on
- `trace-` if USE_CPPTRACE=on

Example:

```bash
mkdir build-deb-mock-tsan && cd build-deb-mock-tsan && \
  ln -s ../data . && \
  cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_TSAN=yes -DUSE_USHER=no .. && \
  make -j32
```

The `data/` symlink is required for tests to find included datasets.

## Testing

Custom test framework. Run from a build directory:

```bash
./bin/larch-test "SPR:.*"           # regex match on test names
./bin/larch-test --range 1,5-10,12  # run specific test IDs
./bin/larch-test +tag slow          # include tests with tag
./bin/larch-test -tag slow          # exclude tests with tag
./bin/larch-test --list             # list matching tests without running
./bin/larch-test nocatch            # let exceptions escape (for debuggers)
```

Never run the complete test suite - it's very slow. Always filter by regex, tags, or ID ranges.

For faster iteration, consider release builds with KEEP_ASSERTS=on when the issue is caught by an assert.

To add a new test:
```cpp
[[maybe_unused]] static const auto registered =
    add_test({[]() { my_test_function(); }, "Test Name", {"optional", "tags"}});
```

## File Formats

**Input:**
- Reference sequences: single-line plain text (optionally gzipped)
- Tree topologies: Newick
- Samples: FASTA, VCF

**Internal/Output:**
- Trees: matOptimize's protobuf format (`.pb.gz`)
- DAGs: Larch's own protobuf format
- JSON: exists but currently unused (may be revived for new features)

## Key Concepts

- **History DAG**: The core data structure. A DAG from which phylogenetic trees can be sampled under various conditions.
- **UA (Universal Ancestor)**: The common root among all trees merged into a DAG.
- **Reference Sequence**: The genome at the UA node.
- **Node/Edge IDs**: Strongly-typed `size_t` used to address nodes and edges.
- **EdgeMutations**: Per-site mutations from parent to child node (stored on edges).
- **CompactGenome**: Per-site mutations from UA to a given node (stored on nodes).
- **SampleID**: String name for leaf nodes identifying the input sequence sample (must be unique).
- **Merging**: Combining trees or fragments to form/expand a DAG.
- **Fragment**: A minimal subset of a tree containing only modified nodes/edges, used for efficient merging.
- **SPR Move**: Subtree Prune and Regraft. Enriches a DAG to hold more trees. Works on a sampled tree with a source and destination node, produces a Fragment.
- **Condensing/Uncondensing**: Technique from matOptimize that groups closely related sibling nodes for performance. Larch condenses before passing to matOptimize and uncondenses the result.

### Scoring and Sampling

Currently only parsimony score is used (1 point per changed base). Complex scoring models are an upcoming major feature.

For sampling strategies and tree scoring, see `include/larch/subtree/*`. CLI arguments for controlling sampling are in `tools/larch-usher.cpp`.

### Merging and Node Identity

Merging is central to how DAGs grow. Key concepts:
- **LeafSet**: The set of leaf descendants, used to identify subtrees.
- **NodeLabel**: Identifies a node by its LeafSet and CompactGenome.
- **EdgeLabel**: Identifies an edge by parent/child NodeLabels and EdgeMutations.

Nodes are considered identical when they have the same NodeLabel. This is how the merge operation detects shared structure between trees/fragments.

See `include/larch/merge/*` and `include/larch/impl/merge/*`.

## Architecture: Storage-View Pattern

The codebase uses an in-house pattern called **storage-view**, analogous to `std::string` vs `std::string_view` but for complex data structures.

- **Storage**: A mixin of a variadic set of Features. Three storage types exist: nodes, edges, and DAG itself.
- **View**: A lightweight handle to a node/edge/DAG. Takes a reference to storage or copies another view (for composability).

This pattern enables layering data and functionality:
- **Overlay**: Store only changes to an existing DAG without modifying the original.
- **Extend**: Add new data and functions to existing structures.

### Combining Features into a DAG

See `include/larch/madag/mutation_annotated_dag.hpp` for how the basic DAG structure combines Features like ReferenceSequence, EdgeMutations, CompactGenomes, and SampleID.

The `LongNameOf<>` helper reduces template names for more readable error messages and stack traces.

## The Optimization Loop

Growing a DAG is iterative:

1. **Sample** a tree from the DAG
2. **Pass to matOptimize**, which floods Larch via `Move_Found_Callback` with candidate SPR moves
3. **Accept/reject** moves individually; accepted moves are merged as Fragments via `BatchingCallback`
4. **Radius exhaustion**: After some moves, matOptimize increases the radius (how far source/destination nodes can be apart) and restarts
5. **Merge whole trees** at the end of each radius
6. **Full iteration complete** when all radii finish; repeat by sampling another tree

## matOptimize Integration

Two ways to cross to/from matOptimize:

- **MAT Conversion** (`include/larch/mat_conversion.hpp`): Deep copy, bidirectional. Creates `MAT::Tree` from MADAG and vice versa.
- **MAT View** (`include/larch/mat_view.hpp`): Lightweight view into `MAT::Tree` that provides a DAG interface, enabling composition with Larch's layered architecture.

Build with USE_MAT_VIEW=off to use only MAT conversion everywhere.

## Parallelism

Both TBB and Taskflow are used. Migration from TBB to Taskflow is ongoing - changes in this direction are welcome. Custom parallel helpers and data structures live in `include/larch/parallel/`.

**Disabling parallelism** (for debugging): Use the `DISABLE_PARALLELISM=yes` cmake option. This:
1. Makes `ParallelForEach()` delegate to `SeqForEach()`
2. Sets `tbb::global_control::max_allowed_parallelism` to 1 at the start of all executables

This disables parallelism throughout Larch's code and in TBB-based components like VCF import pipelines. Note that matOptimize has its own internal TBB usage that respects the global_control setting.

## Debugging

### Classifying bugs

First determine if the bug is:
1. **UB/memory corruption** - run with Address Sanitizer (USE_ASAN=yes)
2. **Threading issue** - disable parallelism first, then use TSan (USE_TSAN=yes)
3. **Logical flaw** - add asserts earlier in the pipeline to narrow down where things go wrong

### Tips

- **Reproduce on small data first**: Debugging on huge trees/DAGs is difficult. Use `test/sample_dag.hpp` which provides a small 10-node tree.
- **Visualization**: Functions `*ToDOT()` in `include/larch/dag_loader.hpp` output Graphviz DOT format. Most useful with small datasets.
- **Stack traces**: Build with USE_CPPTRACE=on for readable backtraces on uncaught exceptions (cheaper than gdb).
- **Thread sanitizer caveats**: TSan with TBB/Taskflow produces many false positives. Use suppressions file and recommended options below.
- **Long-running failures**: When bugs only appear on large datasets after many optimization cycles, adding asserts earlier in the pipeline helps locate the source.

### TSan Configuration

Run TSan builds with these environment options:

```bash
TSAN_OPTIONS="halt_on_error=1 history_size=7 force_seq_cst_atomics=1 second_deadlock_stack=1 suppressions=../tsan.suppressions"
```

The `tsan.suppressions` file in the repo root contains suppressions for known false positives (FreeBSD libthr internals, Taskflow topology cleanup).

### Common Pitfalls

- **Forgetting initialization**: Must call `BuildConnections()` and `ComputeCompactGenomes()` at appropriate times.
- **matOptimize boundary**: matOptimize uses pointers to address nodes (not indices). Helper mappings are needed when crossing this boundary.
- **SPR move confusion**: Easy to mix up pre-move and post-move nodes when working with SPR code.
- **Undocumented matOptimize behavior**: Techniques like `reassign_states` or condensing may have surprising behavior in edge cases.

## Code Layout

- `include/larch/` - header declarations
- `include/larch/impl/` - header definitions (predominantly header-only)
- `include/larch/merge/` - merging logic
- `include/larch/madag/` - mutation-annotated DAG structure
- `include/larch/subtree/` - tree sampling and scoring
- `include/larch/parallel/` - parallel helpers and data structures
- `src/` - non-templated implementation files
- `tools/` - CLI executables (larch-usher, larch-dagutil, larch-dag2dot, bcr-larch)
- `test/` - test suite

## Dependencies

### usher Submodule

`deps/usher` is a git submodule pointing to `git@github.com:matsengrp/usher.git` (tracking `master` branch). Only `deps/usher/src/matOptimize` is used. A mock implementation exists in `include/larch/optimize.hpp` - enable with USE_USHER=no to simplify runtime and rule out matOptimize issues.

**Why not use usher's CMakeLists.txt?** The usher project is not designed to be used as a library. Larch manually compiles a subset of matOptimize sources directly (see `MATOPTIMIZE_SRC` glob in CMakeLists.txt). This avoids pulling in usher's main executable, unused utilities, and conflicting build options.

**Files compiled from matOptimize:**
- `apply_move/*.cpp` - SPR move application
- `Profitable_Moves_Enumerators/*.cpp` - move enumeration
- `Fitch_Sankoff.cpp` - parsimony scoring
- `mutation_annotated_tree*.cpp` - MAT data structure
- `optimize_inner_loop.cpp`, `optimize_tree.cpp` - optimization loop
- `Mutation_Collection.cpp`, `condense.cpp`, `reassign_states.cpp` - tree manipulation
- `detailed_mutations_load.cpp`, `detailed_mutations_store.cpp` - serialization
- `priority_conflict_resolver.cpp`, `check_samples.cpp` - utilities

**Global variables:** matOptimize declares several `extern` globals in `tree_rearrangement_internal.hpp` that are normally defined in `main.cpp`. Since larch doesn't compile `main.cpp`, these definitions are provided in `src/usher_globals.cpp`:
- `interrupted` - atomic flag for interruption handling
- `process_count`, `this_rank` - MPI-related (set to 1/0 for single-process)
- `num_threads` - thread count (initialized from TBB default concurrency)
- `movalbe_src_log` - debug logging file handle
- `use_bound`, `changing_radius` - optimization control flags

**Updating the submodule:** When updating to a new usher commit, check if `tree_rearrangement_internal.hpp` has new `extern` declarations that need definitions added to `usher_globals.cpp`.

## Maintaining This File

This CLAUDE.md is a living document. When working on the codebase, update it with new insights, patterns, or gotchas you discover.

**Areas that could be expanded:**
- SPR move mechanics: more detail on source/destination nodes, how fragments are constructed
- CLI tool usage: what larch-usher, larch-dagutil, larch-dag2dot, bcr-larch actually do
- Performance considerations: memory usage patterns, scaling concerns
- Branch-specific context: the `linearham` branch is for implementing complex scoring models

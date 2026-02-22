# Native MADAG Optimization: Research Notes

Research documentation for reimplementing matOptimize's SPR tree optimization
algorithm to operate directly on larch's native MADAG data structures, without
the MAT conversion roundtrip.

## Documents

1. **[algorithm-overview.md](algorithm-overview.md)** - High-level description of
   matOptimize's optimization loop, radius control, and termination conditions.

2. **[move-enumeration.md](move-enumeration.md)** - Deep dive into the bounded
   search algorithm for finding profitable SPR moves. The most complex and
   performance-critical piece.

3. **[fitch-scoring.md](fitch-scoring.md)** - How parsimony scoring works via
   Fitch sets, major alleles, and boundary1 tracking. Covers both matOptimize's
   incremental approach and larch's existing native Fitch implementation.

4. **[move-application.md](move-application.md)** - How moves are applied to the
   tree: the backward pass (major allele recomputation) and forward pass (state
   propagation).

5. **[conflict-resolution.md](conflict-resolution.md)** - How parallel move
   enumeration resolves conflicts and recycles deferred moves.

6. **[larch-infrastructure.md](larch-infrastructure.md)** - Existing larch
   building blocks that can be reused: SPR overlays, fragments, merge system,
   Fitch backend, and data structures.

7. **[implementation-plan.md](implementation-plan.md)** - Roadmap for
   implementing native optimization, mapping matOptimize concepts to MADAG
   equivalents, and identifying what needs to be built.

8. **[gotchas.md](gotchas.md)** - Pitfalls, edge cases, and non-obvious
   behaviors discovered during research.

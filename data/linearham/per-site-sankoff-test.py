#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "historydag @ git+https://github.com/matsengrp/historydag",
#     "ete3",
#     "numpy",
# ]
# ///
"""Generate ground truth for SankoffScorer validation.

Run: uv run docs/per-site-sankoff-test-plan.py
"""

import numpy as np
import ete3

import historydag.parsimony as dag_parsimony
import historydag.parsimony_utils as parsimony_utils


def create_test_tree():
    """Create test tree with 3 leaves and sequences."""
    tree = ete3.Tree("((L1,L2)n1,L3)root;", format=1)
    seqs = {"L1": "CTAA", "L2": "CACC", "L3": "AGTC"}
    for leaf in tree.iter_leaves():
        leaf.add_feature("sequence", seqs[leaf.name])
    for node in tree.traverse():
        if not node.is_leaf():
            node.add_feature("sequence", "NNNN")
    return tree


def create_matrices():
    """Create 4 site-specific matrices matching C++ test.

    Base indexing: A=0, C=1, G=2, T=3 (ACGT order)
    Matrix indexing: [from_base, to_base]
    """
    matrices = np.zeros((4, 4, 4))

    # Site 1: transitions cheap (1), transversions expensive (2)
    matrices[0] = np.array([
        [0, 2, 1, 2],  # from A
        [2, 0, 2, 1],  # from C
        [1, 2, 0, 2],  # from G
        [2, 1, 2, 0],  # from T
    ])

    # Site 2: uniform
    matrices[1] = np.array([
        [0, 1, 1, 1],
        [1, 0, 1, 1],
        [1, 1, 0, 1],
        [1, 1, 1, 0],
    ])

    # Site 3: A<->G cheap
    matrices[2] = np.array([
        [0, 3, 0.5, 3],
        [3, 0, 3, 3],
        [0.5, 3, 0, 3],
        [3, 3, 3, 0],
    ])

    # Site 4: Asymmetric (A->C=1, C->A=2)
    matrices[3] = np.array([
        [0, 1, 1, 1],  # from A: A->C = 1
        [2, 0, 1, 1],  # from C: C->A = 2 (asymmetric!)
        [1, 1, 0, 1],  # from G
        [1, 1, 1, 0],  # from T
    ])

    return matrices


def compute_per_site_score(tree, site_idx, site_matrix):
    """Compute Sankoff score for a single site by extracting that site's characters.

    Creates a new tree with single-character sequences for just this site.
    """
    # Create a fresh tree with single-character sequences
    t = ete3.Tree("((L1,L2)n1,L3)root;", format=1)
    seqs = {"L1": "CTAA", "L2": "CACC", "L3": "AGTC"}

    for leaf in t.iter_leaves():
        leaf.add_feature("sequence", seqs[leaf.name][site_idx])
    for node in t.traverse():
        if not node.is_leaf():
            node.add_feature("sequence", "N")

    # Use TransitionModel (not Sitewise) for single-site computation
    single_tm = parsimony_utils.TransitionModel(
        bases="ACGT",
        transition_weights=site_matrix
    )

    return dag_parsimony.sankoff_upward(t, seq_len=1, transition_model=single_tm)


def main():
    tree = create_test_tree()
    matrices = create_matrices()

    # IMPORTANT: Must specify bases="ACGT" to match our matrix indexing
    # (historydag defaults to "AGCT" which has different index order)
    tm = parsimony_utils.SitewiseTransitionModel(
        bases="ACGT",
        transition_matrix=matrices
    )

    score = dag_parsimony.sankoff_upward(tree, seq_len=4, transition_model=tm)
    print(f"Total Sankoff score: {score}")
    print(f"Expected total score: 12.0")
    print()

    # Per-site scores (for verification against hand calculations)
    print("Per-site breakdown:")
    site_names = [
        "Site 1 (C,C,A) - ti/tv",
        "Site 2 (T,A,G) - uniform",
        "Site 3 (A,C,T) - A<->G cheap",
        "Site 4 (A,C,C) - asymmetric",
    ]
    expected = [2.0, 2.0, 6.0, 2.0]

    total_from_sites = 0
    for site_idx in range(4):
        site_score = compute_per_site_score(tree, site_idx, matrices[site_idx])
        total_from_sites += site_score
        status = "PASS" if abs(site_score - expected[site_idx]) < 1e-9 else "FAIL"
        print(f"  {site_names[site_idx]}: {site_score} (expected {expected[site_idx]}) {status}")

    print()
    print(f"Sum of per-site scores: {total_from_sites}")
    total_ok = "PASS" if abs(score - 12.0) < 1e-9 else "FAIL"
    print(f"Total score matches expected: {total_ok}")


if __name__ == "__main__":
    main()

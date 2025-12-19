#!/usr/bin/env python3.11

import sys
from pathlib import Path


def reroot_iqtree(tree_string, naive_name):
    """Reroot tree to make naive_name the outgroup."""
    tree_string = tree_string.strip()
    if tree_string.endswith(';'):
        tree_string = tree_string[:-1]

    if not tree_string.startswith('(') or not tree_string.endswith(')'):
        raise ValueError("Invalid tree format")

    # Remove outer parens to get inner content
    inner = tree_string[1:-1]

    # Split by top-level commas (depth 0)
    children = []
    depth = 0
    start = 0
    for i, c in enumerate(inner):
        if c == '(':
            depth += 1
        elif c == ')':
            depth -= 1
        elif c == ',' and depth == 0:
            children.append(inner[start:i])
            start = i + 1
    children.append(inner[start:])  # Last child

    # Find the naive sequence child by name
    naive_child = None
    naive_idx = None
    for i, child in enumerate(children):
        # Check if this child is the naive sequence (simple leaf: "name:branch_length")
        if child.startswith(naive_name + ':') or child == naive_name:
            naive_child = child
            naive_idx = i
            break

    if naive_child is None:
        raise ValueError(f"Could not find naive sequence '{naive_name}' at top level")

    # Extract branch length from naive
    if ':' in naive_child:
        _, branch_length = naive_child.split(':', 1)
    else:
        branch_length = "0.0"

    # Remove naive from children and reconstruct
    # Make naive the label of the root node (internal node below UA)
    remaining_children = children[:naive_idx] + children[naive_idx + 1:]
    rest_of_tree = ','.join(remaining_children)

    # The naive sequence becomes the root node label, not a leaf
    return f"({rest_of_tree}){naive_name}:{branch_length};"


def get_sequence_names(fasta_path):
    """Extract sequence names from a FASTA file."""
    names = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                names.append(line[1:].strip())
    return names


def main():
    naive_fa = Path("naive-seqs.fa")

    if not naive_fa.exists():
        print(f"Error: {naive_fa} not found", file=sys.stderr)
        sys.exit(1)

    seq_names = get_sequence_names(naive_fa)

    for name in seq_names:
        treefile = Path(f"{name}.treefile")

        if not treefile.exists():
            print(f"Warning: {treefile} not found, skipping", file=sys.stderr)
            continue

        tree_string = treefile.read_text().strip()
        rerooted = reroot_iqtree(tree_string, name)
        output_file = Path(f"{name}-rerooted.treefile")
        output_file.write_text(rerooted + "\n")
        print(f"Re-rooted {treefile} -> {output_file}")


if __name__ == "__main__":
    main()

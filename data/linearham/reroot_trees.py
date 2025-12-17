#!/usr/bin/env python3.11

import sys
from pathlib import Path


def reroot_iqtree(tree_string):
    first_comma = tree_string.find(',')
    root_info = tree_string[1:first_comma]
    root_name, branch_length = root_info.split(':')
    rest_of_tree = tree_string[first_comma+1:-1]
    new_tree_string = "((" + rest_of_tree + ":" + branch_length + ")" + root_name + ";"

    return new_tree_string


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
        rerooted = reroot_iqtree(tree_string)
        output_file = Path(f"{name}-rerooted.treefile")
        output_file.write_text(rerooted + "\n")
        print(f"Re-rooted {treefile} -> {output_file}")


if __name__ == "__main__":
    main()

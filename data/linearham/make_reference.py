#!/usr/bin/env python3.11

from collections import Counter
from pathlib import Path


def read_sequences(fasta_path):
    """Read sequences from a FASTA file."""
    sequences = []
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)

        if current_seq:
            sequences.append(''.join(current_seq))

    return sequences


def majority_rule_consensus(sequences):
    """Create a consensus sequence using majority rule at each position."""
    if not sequences:
        return ""

    max_len = max(len(seq) for seq in sequences)
    consensus = []

    for i in range(max_len):
        bases = [seq[i] for seq in sequences if i < len(seq)]
        counter = Counter(bases)
        most_common = counter.most_common(1)[0][0]
        consensus.append(most_common)

    return ''.join(consensus)


def main():
    naive_fa = Path("naive-seqs.fa")
    output_file = Path("reference_sequence.txt")

    sequences = read_sequences(naive_fa)
    consensus = majority_rule_consensus(sequences)
    output_file.write_text(consensus + "\n")

    print(f"Created {output_file} ({len(consensus)} bases)")


if __name__ == "__main__":
    main()

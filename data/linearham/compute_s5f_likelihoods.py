#!/usr/bin/env python3.11
"""
Compute per-edge log-likelihoods using the S5F model for a test tree.

This script generates ground truth values that can be compared with the C++
implementation in larch.
"""

import sys
import torch
import yaml

# Reference sequence from linearham data (381 bases)
REFERENCE_SEQ = (
    "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTC"
    "TGGATTCACCGTCAGTAGCAACTACATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAG"
    "TTATTTATAGCGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCAGATTCACCATCTCCAGAGACAATTCC"
    "AAGAACACGCTGTATCTTCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGGCAC"
    "AACACACGGGTATAGCAGTGAAGGCATGACTTCAAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCG"
    "TCTCCTCAG"
)


def generate_kmers(k):
    """Generate all possible k-mers."""
    bases = ['A', 'C', 'G', 'T']
    if k == 0:
        return ['']
    smaller = generate_kmers(k - 1)
    return [b + s for b in bases for s in smaller]


def kmer_to_index(kmer):
    """Convert a k-mer to its index."""
    base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    index = 0
    for b in kmer:
        if b not in base_map:
            return -1  # N or invalid
        index = index * 4 + base_map[b]
    return index


def encode_sequence(seq, kmer_length, site_count):
    """Encode sequence as kmer indices."""
    half_k = kmer_length // 2
    padded = 'N' * half_k + seq + 'N' * half_k
    invalid_idx = 4 ** kmer_length  # Last index for N-containing kmers

    encoded = []
    for i in range(len(seq)):
        kmer = padded[i:i + kmer_length]
        idx = kmer_to_index(kmer)
        if idx < 0:
            idx = invalid_idx  # Use last index for invalid kmers
        encoded.append(idx)

    # Pad or truncate to site_count
    while len(encoded) < site_count:
        encoded.append(0)
    encoded = encoded[:site_count]

    return torch.tensor(encoded, dtype=torch.long)


def create_wt_modifier(seq, site_count):
    """Create wildtype base modifier tensor."""
    base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    modifier = torch.zeros(site_count, 4)

    for i, base in enumerate(seq[:site_count]):
        if base in base_map:
            modifier[i, base_map[base]] = float('-inf')

    return modifier


def create_mask(seq, site_count):
    """Create mask tensor for valid positions."""
    mask = torch.ones(site_count, dtype=torch.float32)
    for i, base in enumerate(seq[:site_count]):
        if base == 'N':
            mask[i] = 0.0
    # Zero out positions beyond sequence length
    for i in range(len(seq), site_count):
        mask[i] = 0.0
    return mask


class RSFivemerModel(torch.nn.Module):
    """Pure PyTorch implementation of RSFivemerModel."""

    def __init__(self, kmer_count):
        super().__init__()
        self.r_kmer_embedding = torch.nn.Embedding(kmer_count, 1)
        self.s_kmer_embedding = torch.nn.Embedding(kmer_count, 4)

    def forward(self, encoded_parents, masks, wt_base_modifier):
        # R component: get log rates from embedding and exponentiate
        log_kmer_rates = self.r_kmer_embedding(encoded_parents).squeeze(-1)
        rates = torch.exp(log_kmer_rates * masks)

        # S component: get CSP logits from embedding
        csp_logits = self.s_kmer_embedding(encoded_parents)

        # When we have an N, set all the CSP logits to 0
        csp_logits = csp_logits * masks.unsqueeze(-1)

        # Apply wt_base_modifier to zero out wild-type base predictions
        csp_logits = csp_logits + wt_base_modifier

        return rates, csp_logits


def load_model(pth_path, yml_path):
    """Load the S5F model."""
    with open(yml_path, 'r') as f:
        config = yaml.safe_load(f)

    kmer_length = config['encoder_parameters']['kmer_length']
    site_count = config['encoder_parameters']['site_count']
    # Extra embedding for N-containing kmers (index 4^k)
    kmer_count = 4 ** kmer_length + 1

    model = RSFivemerModel(kmer_count)

    # Load state dict
    state_dict = torch.load(pth_path, weights_only=True)
    model.load_state_dict(state_dict)
    model.eval()

    return model, kmer_length, site_count


def mutate_sequence(seq: str, pos: int, new_base: str) -> str:
    """Create a new sequence with a single mutation."""
    return seq[:pos] + new_base + seq[pos + 1:]


def compute_log_likelihood(model, kmer_length, site_count,
                          parent_seq: str, child_seq: str) -> float:
    """
    Compute the log-likelihood for a parent-child edge using the Poisson
    context model.

    This matches the formula in netam:
    log_lik = sum(log(lambda_j)) + n * log(t_hat) - n
    where:
    - lambda_j = rate_j * csp_j (mutation probability at position j)
    - n = number of mutations
    - t_hat = n / sum(rates) (estimated evolution time)
    """
    with torch.no_grad():
        # Encode parent sequence
        encoded = encode_sequence(parent_seq, kmer_length, site_count).unsqueeze(0)
        wt_mod = create_wt_modifier(parent_seq, site_count).unsqueeze(0)
        mask = create_mask(parent_seq, site_count).unsqueeze(0)

        # Run forward pass
        rates, csp_logits = model(encoded, mask, wt_mod)

        # Apply softmax to get CSP probabilities
        csp = torch.softmax(csp_logits, dim=-1)

        # Encode base indices (A=0, C=1, G=2, T=3)
        base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
        parent_indices = [base_to_idx.get(b, 4) for b in parent_seq]
        child_indices = [base_to_idx.get(b, 4) for b in child_seq]

        # Trim to site_count
        seq_len = min(len(parent_seq), site_count)
        parent_indices = parent_indices[:seq_len]
        child_indices = child_indices[:seq_len]

        # Find mutation positions
        n_mutations = 0
        log_lambda_sum = 0.0
        rates_squeezed = rates.squeeze(0)
        csp_squeezed = csp.squeeze(0)

        for pos in range(seq_len):
            p_idx = parent_indices[pos]
            c_idx = child_indices[pos]
            if p_idx != c_idx and p_idx < 4 and c_idx < 4:
                n_mutations += 1
                rate_j = rates_squeezed[pos].item()
                csp_j = csp_squeezed[pos, c_idx].item()
                lambda_j = rate_j * csp_j
                if lambda_j > 0:
                    log_lambda_sum += torch.log(torch.tensor(lambda_j)).item()
                else:
                    log_lambda_sum += float('-inf')

        if n_mutations == 0:
            return 0.0

        # Compute t_hat = n / sum(rates)
        mask_squeezed = mask.squeeze(0)
        sum_rates = (rates_squeezed * mask_squeezed).sum().item()
        if sum_rates <= 0:
            return float('-inf')

        t_hat = n_mutations / sum_rates

        # Final log-likelihood: sum(log(lambda_j)) + n * log(t_hat) - n
        log_lik = log_lambda_sum + n_mutations * torch.log(torch.tensor(t_hat)).item() - n_mutations

        return log_lik


def main():
    # Load S5F model
    pth_path = "/home/ogi/matsen/larch-s5f-likelihood/data/linearham/s5f.pth"
    yml_path = "/home/ogi/matsen/larch-s5f-likelihood/data/linearham/s5f.yml"
    model, kmer_length, site_count = load_model(pth_path, yml_path)

    # Print sequence info for debugging
    print(f"// Reference sequence length: {len(REFERENCE_SEQ)}")

    # Define test edges with actual mutations based on what's in the sequence
    # Positions and their actual bases (0-indexed):
    # Position 10: T
    # Position 50: C
    # Position 100: T
    # Position 150: A
    # Position 200: C
    test_edges = []

    # Verify bases are correct
    print(f"// Position 10: {REFERENCE_SEQ[10]}")
    print(f"// Position 50: {REFERENCE_SEQ[50]}")
    print(f"// Position 100: {REFERENCE_SEQ[100]}")
    print(f"// Position 150: {REFERENCE_SEQ[150]}")
    print(f"// Position 200: {REFERENCE_SEQ[200]}")
    print()

    # Edge 0: root (ref) -> mutation at position 10: T->A
    child_1 = mutate_sequence(REFERENCE_SEQ, 10, 'A')
    test_edges.append(("edge_1_T_to_A_pos10", REFERENCE_SEQ, child_1, 10, 'T', 'A'))

    # Edge 1: root (ref) -> mutation at position 50: C->T
    child_2 = mutate_sequence(REFERENCE_SEQ, 50, 'T')
    test_edges.append(("edge_2_C_to_T_pos50", REFERENCE_SEQ, child_2, 50, 'C', 'T'))

    # Edge 2: root (ref) -> mutation at position 100: T->G
    child_3 = mutate_sequence(REFERENCE_SEQ, 100, 'G')
    test_edges.append(("edge_3_T_to_G_pos100", REFERENCE_SEQ, child_3, 100, 'T', 'G'))

    # Edge 3: root (ref) -> mutation at position 150: A->C
    child_4 = mutate_sequence(REFERENCE_SEQ, 150, 'C')
    test_edges.append(("edge_4_A_to_C_pos150", REFERENCE_SEQ, child_4, 150, 'A', 'C'))

    # Edge 4: root (ref) -> mutation at position 200: C->A
    child_5 = mutate_sequence(REFERENCE_SEQ, 200, 'A')
    test_edges.append(("edge_5_C_to_A_pos200", REFERENCE_SEQ, child_5, 200, 'C', 'A'))

    # Edge 5: identical sequences (no mutations)
    test_edges.append(("edge_6_identical", REFERENCE_SEQ, REFERENCE_SEQ, -1, '', ''))

    # Edge 6: multiple mutations at positions 10, 50, 100
    # Positions: 10=T, 50=C, 100=T
    multi_mut = mutate_sequence(REFERENCE_SEQ, 10, 'G')  # T->G
    multi_mut = mutate_sequence(multi_mut, 50, 'A')       # C->A
    multi_mut = mutate_sequence(multi_mut, 100, 'C')      # T->C
    test_edges.append(("edge_7_multi_mutation", REFERENCE_SEQ, multi_mut, -1, '', ''))

    print("// Ground truth log-likelihoods computed with S5F model from Python")
    print("// These values should match the C++ implementation")
    print()
    print("struct S5FGroundTruth {")
    print("  const char* name;")
    print("  double log_likelihood;")
    print("};")
    print()
    print("const S5FGroundTruth kS5FGroundTruth[] = {")

    for entry in test_edges:
        name, parent, child = entry[0], entry[1], entry[2]
        ll = compute_log_likelihood(model, kmer_length, site_count, parent, child)
        print(f'    {{"{name}", {ll:.15g}}},')

    print("};")
    print()

    # Also print for easy verification
    print("// Summary:")
    for entry in test_edges:
        name, parent, child = entry[0], entry[1], entry[2]
        ll = compute_log_likelihood(model, kmer_length, site_count, parent, child)
        if len(entry) > 3 and entry[3] >= 0:
            pos, old_base, new_base = entry[3], entry[4], entry[5]
            print(f"// {name} (pos {pos}: {old_base}->{new_base}): {ll:.15g}")
        else:
            print(f"// {name}: {ll:.15g}")


if __name__ == "__main__":
    main()

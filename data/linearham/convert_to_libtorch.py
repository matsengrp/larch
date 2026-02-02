#!/usr/bin/env python3.11
"""
Convert a PyTorch model saved with OrderedDict state_dict to a format
compatible with libtorch's torch::pickle_load (plain dict).

Usage:
    python3.11 convert_to_libtorch.py input.pth output-libtorch.pth

If no output path is given, appends '-libtorch' before the extension:
    python3.11 convert_to_libtorch.py s5f.pth
    # Creates s5f-libtorch.pth
"""

import sys
import torch
from pathlib import Path


def convert_to_libtorch(input_path: str, output_path: str = None) -> str:
    """
    Convert a PyTorch state_dict to libtorch-compatible format.

    Args:
        input_path: Path to the input .pth file
        output_path: Path for the output file (optional)

    Returns:
        Path to the output file
    """
    input_path = Path(input_path)

    if output_path is None:
        output_path = input_path.with_stem(input_path.stem + "-libtorch")
    else:
        output_path = Path(output_path)

    # Load the state dict
    state = torch.load(input_path, weights_only=True)

    # Convert OrderedDict to plain dict if necessary
    if hasattr(state, 'keys'):
        state_dict = dict(state)
    else:
        raise ValueError(f"Expected a state dict, got {type(state)}")

    # Print info about the conversion
    print(f"Input: {input_path}")
    print(f"Output: {output_path}")
    print(f"Original type: {type(state).__name__}")
    print(f"Converted type: {type(state_dict).__name__}")
    print(f"Keys:")
    for key in state_dict.keys():
        tensor = state_dict[key]
        print(f"  {key}: {tuple(tensor.shape)}")

    # Save in libtorch-compatible format
    torch.save(state_dict, output_path)

    print(f"\nSaved to {output_path}")
    return str(output_path)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else None

    convert_to_libtorch(input_path, output_path)


if __name__ == "__main__":
    main()

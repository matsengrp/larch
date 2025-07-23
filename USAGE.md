# Larch Tools Usage Guide

This document provides usage instructions for the executables in the `tools/` directory.

## larch-usher

Tool for exploring tree space of DAG/tree through SPR moves.

### Usage
```shell
larch-usher -i,--input FILE -o,--output FILE [-c,--count INT]
```

### Options
- `-i,--input FILE` - Path to input DAG/Tree file (REQUIRED)
- `-o,--output FILE` - Path to output DAG file (REQUIRED)
- `-r,--MAT-refseq-file FILE` - Path to json reference sequence file (REQUIRED if input file is a MAT protobuf)
- `-v,--VCF-input-file FILE` - Path to VCF file, containing ambiguous leaf sequence data
- `-c,--count INT` - Number of iterations (default: 1)
- `--inter-save INT` - Saves a new intermediate DAG file once every given number of iterations (default: no intermediate DAG files saved)
- `-s,--switch-subtrees INT` - Switch to optimizing subtrees after the specified number of iterations (default: never)
- `--min-subtree-clade-size INT` - The minimum number of leaves in a subtree sampled for optimization (default: 100, ignored without option `-s`)
- `--max-subtree-clade-size INT` - The maximum number of leaves in a subtree sampled for optimization (default: 1000, ignored without option `-s`)
- `--move-coeff-nodes INT` - New node coefficient for scoring moves (default: 1)
- `--move-coeff-pscore INT` - Parsimony score coefficient for scoring moves (default: 1)
- `--sample-any-tree` - Sample any tree for optimization, rather than requiring the sampled tree to maximize parsimony
- `--sample-uniformly` - Use a uniform distribution to sample trees for optimization, rather than a natural distribution
- `--sample-method ENUM` - Select tree sampling method for optimization (default: max parsimony)
  - Options: `parsimony`, `random`, `rf-minsum`, `rf-maxsum`
- `--callback-option ENUM` - Callback configuration choice (default: merge all profitable moves)
  - Options: `best-move`, `best-move-fixed-tree`, `best-move-treebased`, `all-moves`
- `--trim` - Trim optimized DAG after final iteration
- `--keep-fragment-uncollapsed` - Keep empty fragment edges, rather than collapsing them
- `--input-format ENUM` - Specify input file format (default: inferred)
  - Options: `dagbin`, `dag-pb`, `tree-pb`, `dag-json`
- `--output-format ENUM` - Specify output file format (default: inferred)
  - Options: `dagbin`, `dag-pb`
- `--seed INT` - Set seed for random number generation (default: random)
- `--thread INT` - Set number of cpu threads (default: max allowed by system)
- `-T,--max-time` - Exit after fixed runtime (in minutes)
- `-S,--autodetect-stoptime` - Set program to exit after parsimony improvement plateaus
- `-h,--help` - Show help message
- `--version` - Show version information

### File Formats

The tools support several file formats:

- **dagbin** - Binary DAG format
- **dag-pb** - Protobuf DAG format
- **tree-pb** - Protobuf tree format
- **dag-json** - JSON DAG format
- **txt** - Plain text format for reference sequence file

When format is not specified, the tools will attempt to infer the format from the file extension.
Note: larch-usher writes output in either DAG protobuf format or else in Binary DAG format.

### Sample Methods

The `--sample-method` option controls how trees are sampled from the DAG for optimization:

- **parsimony** (default) - Sample trees that maximize parsimony
- **random** - Sample any tree randomly
- **rf-minsum** - Sample trees that minimize the sum of Robinson-Foulds distances
- **rf-maxsum** - Sample trees that maximize the sum of Robinson-Foulds distances

Note: The `random` and `parsimony` methods support the additional `--sample-uniformly` option to use uniform distribution instead of natural distribution.

### Callback Options

The `--callback-option` controls the move selection strategy:

- **best-move** (default) - Merge all profitable moves
- **best-move-fixed-tree** - Find best moves on a fixed tree
- **best-move-treebased** - Tree-based move selection
- **all-moves** - Accept all moves regardless of score

### Examples
```shell
# Basic optimization with 10 iterations
larch-usher -i input.dagbin -o optimized.dagbin -c 10
# MB Alternative command with specific data:
larch-usher -i data/testcase/tree_0.pb.gz -o optimized.dagbin -c 10

# Optimize MAT protobuf with reference sequence
larch-usher -i tree.pb -r refseq.json -o optimized.dagbin -c 5
# MB Alternative command with specific data:
larch-usher -i data/seedtree/seedtree.pb.gz -r data/seedtree/refseq.txt.gz -o optimized.dagbin -c 5

# Use VCF data for ambiguous sequences
larch-usher -i input.dagbin -v sequences.vcf -o output.dagbin
# MB Alternative command with specific data:
larch-usher -i data/test_ambiguous_vcf/amb_mat.pb -r data/test_ambiguous_vcf/sample_reference_sequence.fasta -v data/test_ambiguous_vcf/amb.vcf -o output.dagbin

# Optimize with subtree switching after 20 iterations
larch-usher -i input.dagbin -o output.dagbin -c 50 -s 20
# MB Alternative command with specific data:
larch-usher -i data/testcase/full_dag.pb -o output.dagbin -c 50 -s 20

# Optimize with subtree switching after 20 iterations and setting maximum number of leaves to 20 for sampled subtrees
larch-usher -i input.dagbin -o output.dagbin -c 50 -s 20 --max-subtree-clade-size 20
# MB Alternative command with specific data:
larch-usher -i data/testcase/full_dag.pb -o output.dagbin -c 50 -s 20 --max-subtree-clade-size 20

# Save intermediate results every 5 iterations
larch-usher -i input.dagbin -o final.dagbin -c 20 --inter-save 5
# MB Alternative command with specific data:
larch-usher -i data/testcase/full_dag.pb -o final.dagbin -c 20 --inter-save 5

# Use specific sampling method and callback
larch-usher -i input.dagbin -o output.dagbin --sample-method rf-minsum --callback-option best-move
# MB Alternative command with specific data:
larch-usher -i data/testcase/full_dag.pb -o output.dagbin --sample-method rf-minsum --callback-option best-move

# Set time limit and use auto-stopping
larch-usher -i input.dagbin -o output.dagbin -T 60 -S
# MB Alternative command with specific data:
larch-usher -i data/testcase/full_dag.pb -o output.dagbin -T 60 -S

# Custom move scoring coefficients
larch-usher -i input.dagbin -o output.dagbin --move-coeff-nodes 2 --move-coeff-pscore 3
# MB Alternative command with specific data:
larch-usher -i data/testcase/full_dag.pb -o output.dagbin --move-coeff-nodes 2 --move-coeff-pscore 3

# Final trimming and custom thread count
larch-usher -i input.dagbin -o output.dagbin -c 10 --trim --thread 8
# MB Alternative command with specific data:
larch-usher -i data/testcase/full_dag.pb -o output.dagbin -c 10 --trim --thread 8

# Complete workflow: merge and optimize trees
larch-dagutil -i tree1.dagbin -i tree2.dagbin -i tree3.dagbin -o merged.dagbin
larch-usher -i merged.dagbin -o optimized.dagbin -c 20
larch-dagutil -i optimized.dagbin -t -o final.dagbin
# MB Alternative command with specific data:
larch-dagutil -i data/testcase/tree_1.pb.gz -i data/testcase/tree_2.pb.gz -i data/testcase/tree_3.pb.gz -o merged.dagbin
larch-usher -i merged.dagbin -o optimized.dagbin -c 20
larch-dagutil -i optimized.dagbin -t -o final.dagbin
```

## larch-dagutil

General utility for manipulating (e.g. combining, pruning) or inspecting DAGs/trees.

### Usage
```shell
larch-dagutil [-r,--refseq FILE] -i,--input FILE1 [-i,--input FILE2] ... [-o,--output FILE]
```

### Options
- `-i,--input FILE ` - Paths to input DAG/Tree files (REQUIRED)
- `-o,--output FILE` - Path to output DAG file (default: does not save result DAG)
- `-r,--MAT-refseq-file FILE` - Path to json reference sequence file (REQUIRED if input file is a MAT protobuf)
- `-t,--trim` - Trim output (default: best parsimony)
- `--rf FILE` - Trim output to minimize RF distance to provided DAG file
- `-s,--sample` - Sample a single tree from DAG
- `--dag-info` - Print DAG info (parsimony scores, sum RF distances)
- `--parsimony` - Print all DAG parsimony scores
- `--sum-rf-distance` - Print all DAG sum RF distances
- `--input-format ENUM [...]` - Specify input file formats (default: inferred)
  - Options: `dagbin`, `dag-pb`, `tree-pb`, `dag-json` (see larch-usher for file format details)
- `--output-format ENUM` - Specify output file format (default: inferred)
  - Options: `dagbin`, `dag-pb`
- `--rf-format ENUM` - Specify RF file format (default: inferred)
  - Options: `dagbin`, `dag-pb`, `tree-pb`, `dag-json`
- `-h,--help` - Show help message
- `--version` - Show version information

### Examples
```shell
# Merge multiple DAG files
larch-dagutil -i dag1.dagbin dag2.dagbin -o merged.dagbin

# Merge trees with reference sequence
larch-dagutil -r refseq.json -i tree1.pb -i tree2.pb -o merged_dag.pb
# MB Alternative command with specific data:
larch-dagutil -r data/seedtree/refseq.txt.gz -i data/seedtree/seedtree.pb.gz -i data/seedtree/seedtree_optimized_MAT.pb -o merged_dag.pb

# Trim merged DAG to best parsimony
larch-dagutil -i dag1.dagbin -i dag2.dagbin -t -o trimmed.dagbin
# MB Alternative command with specific data:
larch-dagutil -i data/testcase/tree_0.pb.gz -i data/testcase/tree_1.pb.gz -t -o trimmed.dagbin

# Sample a single tree from DAG
larch-dagutil -i input.dagbin -s -o sampled_tree.dagbin
# MB Alternative command with specific data:
larch-dagutil -i data/testcase/full_dag.pb -s -o sampled_tree.dagbin

# Print DAG information without saving
larch-dagutil -i input.dagbin --dag-info
# MB Alternative command with specific data:
larch-dagutil -i data/testcase/full_dag.pb --dag-info

# Trim to minimize RF distance to reference DAG
larch-dagutil -i input.dagbin --rf reference.dagbin -t -o trimmed.dagbin
# MB Alternative command with specific data:
larch-dagutil -i data/testcase/full_dag.pb --rf data/testcase/tree_0.pb.gz -o trimmed.dagbin

# Convert MAT to DAG
larch-dagutil -r refseq.json -i tree.pb -o dag.dagbin
# MB Alternative command with specific data:
larch-dagutil -r data/seedtree/refseq.txt.gz -i data/seedtree/seedtree.pb.gz -o dag.dagbin

# Complete workflow: merge multiple trees and trim to best parsimony
larch-dagutil -i tree1.dagbin -i tree2.dagbin -i tree3.dagbin -o merged.dagbin
larch-dagutil -i merged.dagbin -t -o final.dagbin
# MB Alternative command with specific data:
larch-dagutil -i data/testcase/tree_0.pb.gz -i data/testcase/tree_1.pb.gz -i data/testcase/tree_2.pb.gz -o merged.dagbin
larch-dagutil -i merged.dagbin -t -o final.dagbin
```

## larch-dag2dot

Converts DAG/tree files to DOT format for visualization.

### Usage
```shell
larch-dag2dot -i,--input FILE [-o,--output FILE]
```

### Options
- `-i,--input FILE` - Path to input DAG/Tree file (REQUIRED)
- `-o,--output FILE` - Path to output DOT file (default: DOT written to stdout)
- `--input-format ENUM` - Specify input file format (default: inferred)
  - Options: `dagbin`, `dag-pb`, `tree-pb`, `dag-json` (see larch-usher for file format details)
- `--dag/--tree` - Specify whether protobuf input is a DAG or Tree
- `-h,--help` - Show help message
- `--version` - Show version information

### Examples
```shell
# Convert a DAG binary file to DOT format and write to stdout
larch-dag2dot -i input.dagbin
# MB Alternative command with specific data:
larch-dag2dot -i data/testcase/full_dag.pb

# Convert a tree protobuf to DOT and save to file
larch-dag2dot -i tree.pb --tree -o output.dot
# MB Alternative command with specific data:
larch-dag2dot -i data/seedtree/seedtree.pb.gz --tree -o output.dot

# Specify input format explicitly
larch-dag2dot -i input.dag --input-format dag-json -o graph.dot
# MB Alternative command with specific data:
larch-dag2dot -i data/test_5_trees/full_dag.json.gz --input-format dag-json -o graph.dot

# Visualize a DAG as PNG image
larch-dag2dot -i mydag.dagbin | dot -Tpng -o mydag.png
# MB Alternative command with specific data:
larch-dag2dot -i data/test_5_trees/full_dag.pb.gz | dot -Tpng -o mydag.png
```

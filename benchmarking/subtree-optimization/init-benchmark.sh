#!/bin/bash
set -eu

./run-option-benchmark.sh 20D_intermediate ../../data/20D_from_fasta/1final-tree-1.nh1.pb.gz ../../data/20D_from_fasta/refseq.txt.gz 1000000

# ./run-option-benchmark.sh startmat ../../data/startmat/startmat_no_ancestral.pb.gz ../../data/startmat/refseq.txt.gz 1000000

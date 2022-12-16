#!/bin/bash
set -eu

./run-option-benchmark.sh 20D/log ../../data/20D_from_fasta/1final-tree-1.nh1.pb.gz ../../data/20D_from_fasta/refseq.txt.gz 5000

./run-option-benchmark.sh startmat/log ../../data/startmat/startmat_no_ancestral.pb.gz ../../data/startmat/refseq.txt.gz 5000

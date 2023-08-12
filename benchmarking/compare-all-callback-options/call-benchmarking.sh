#!/bin/bash
set -eu
packagepath="benchmark-callbacks"
fp="../../data/seedtree/seedtree_optimized_MAT.pb"
rsfp="../../data/seedtree/refseq.txt.gz"
numiter=50
numreps=50

./benchmark.sh $packagepath $fp $rsfp $numiter $numreps

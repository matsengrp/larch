#!/bin/bash
set -eu
packagepath="benchmark-callbacks"
fp="../../data/seedtree/seedtree.pb.gz"
rsfp="../../data/seedtree/refseq.txt.gz"
numiter=50
numreps=50

./benchmark.sh $packagepath $fp $rsfp $numiter $numreps

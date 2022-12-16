#!/bin/bash
set -eu
# eval "$(conda shell.bash hook)"
# conda activate usher
logprefix=$1
inmat=$2
refseq=$3
iterations=$4
numruns=50


for option1 in {'','--sample-best-tree',}; do
    for option2 in {'','--sample-uniform',}; do
        logfilepath=$logprefix/option_$option1$option2
        echo $logfilepath
        mkdir -p $logfilepath
        for jobnum in $(seq $numruns); do
            larch_usher_options="-i $inmat -o opt_dag.pb -c $iterations -r $refseq -l $logfilepath/log_$jobnum $option1 $option2"
            sbatch -c 4 -J log$jobnum -o log$jobnum.log ./run-larch-usher.sh $larch_usher_options
        done
    done
done



#!/bin/bash
set -eu
# eval "$(conda shell.bash hook)"
# conda activate usher
logprefix=$1
inmat=$2
refseq=$3
iterations=$4
numruns=50


for option1 in "--sample-best-tree" "--sample-best-tree -s 0" "--sample-best-tree -s 2" "--sample-best-tree -s 10"; do
    optionid=${option1//[[:blank:]]/}
    logfilepath=$logprefix/option_$optionid
    echo $logfilepath
    mkdir -p $logfilepath
    for jobnum in $(seq $numruns); do
        larch_usher_options="-i $inmat -o opt_dag.pb -c $iterations -r $refseq -l $logfilepath/log_$jobnum $option1"
        sbatch -c 4 -J rn$jobnum -o $logfilepath/log$jobnum.log ../run-larch-usher.sh $larch_usher_options
    done
done

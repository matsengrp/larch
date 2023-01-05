#!/bin/bash
set -eu
# eval "$(conda shell.bash hook)"
# conda activate usher
logprefix=$1
inmat=$2
refseq=$3
iterations=$4
numruns=20


option2="--min-subtree-clade-size 100 --max-subtree-clade-size 1000 --uniform-subtree-root"
option3="--min-subtree-clade-size 2 --max-subtree-clade-size 70000 --uniform-subtree-root"

for option1 in "" "-s 0 $option3" "-s 0 $option2" "-s 10 $option3" "-s 10 $option2"; do
    optionid=${option1//[[:blank:]]/}
    logfilepath=$logprefix/option_$optionid
    echo $logfilepath
    mkdir -p $logfilepath
    for jobnum in $(seq $numruns); do
        larch_usher_options="-i $inmat -o opt_dag.pb -c $iterations -r $refseq -l $logfilepath/log_$jobnum $option1"
        sbatch -c 2 -J rn$jobnum -o $logfilepath/log$jobnum.log ../run-larch-usher.sh $larch_usher_options
    done
done

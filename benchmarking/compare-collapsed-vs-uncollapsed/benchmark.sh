#!/bin/bash
set -eu
# eval "$(conda shell.bash hook)"
# conda activate usher
logprefix=$1
inmat=$2
refseq=$3
iterations=$4
numruns=$5

option1names=('collapse-fragment' 'uncollapse-fragment')
option1args=('--callback-option best-moves' '--callback-option best-moves --keep-fragment-uncollapsed')
for i in "${!option1names[@]}"; do
    option1="${option1names[i]}"
    option1arg="${option1args[i]}"
    logfilepath=$logprefix/option_$option1
    echo $logfilepath
    mkdir -p $logfilepath
    for jobnum in $(seq $numruns); do
        larch_usher_options="-i $inmat -o opt_dag.pb -c $iterations -r $refseq -l $logfilepath/log_$jobnum $option1arg"
        sbatch -c 4 -J callbackoption$jobnum -o $jobnum.log ../run-larch-usher.sh $larch_usher_options
    done
done

pythonparams=("$inmat" "$iterations" "$numruns" "$logprefix" "+" "${option1names[@]}")
printf "%s\n" "${pythonparams[@]}" > python_params.txt


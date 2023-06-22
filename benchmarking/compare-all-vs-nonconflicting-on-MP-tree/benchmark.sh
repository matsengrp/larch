#!/bin/bash
set -eu
# eval "$(conda shell.bash hook)"
# conda activate usher
logprefix=$1
inmat=$2
refseq=$3
iterations=$4
numruns=$5

optionlist=('best-moves' 'best-moves-treebased')
for option in ${optionlist[@]}; do
    logfilepath=$logprefix/option_$option
    echo $logfilepath
    mkdir -p $logfilepath
    for jobnum in $(seq $numruns); do
        larch_usher_options="-i $inmat -o opt_dag.pb -c $iterations -r $refseq -l $logfilepath/log_$jobnum --callback-option $option"
        sbatch -c 4 -J callbackoption$jobnum -o $jobnum.log ../run-larch-usher.sh $larch_usher_options
    done
done

pythonparams=("$inmat" "$iterations" "$numruns" "$logprefix" "+" "${optionlist[@]}")
printf "%s\n" "${pythonparams[@]}" > python_params.txt


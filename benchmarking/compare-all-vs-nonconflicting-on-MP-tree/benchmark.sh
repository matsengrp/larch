#!/bin/bash
set -eu
# eval "$(conda shell.bash hook)"
# conda activate usher
logprefix=$1
inmat=$2
refseq=$3
iterations=$4
numruns=$5

optionargs=('best-moves' 'best-moves-treebased' 'best-moves --keep-fragment-uncollapsed')
optionnames=('all-profitable-moves' 'nonconflicting-profitable-moves' 'all-profitable-moves-uncollapsed')
for optioni in ${!optionargs[@]}; do
    option="${optionargs[optioni]}"
    optionname="${optionnames[optioni]}"
    logfilepath=$logprefix/option_$optionname
    echo $logfilepath
    mkdir -p $logfilepath
    for jobnum in $(seq $numruns); do
        larch_usher_options="-i $inmat -o opt_dag.pb -c $iterations -r $refseq -l $logfilepath/log_$jobnum --callback-option $option"
        sbatch -c 4 -J callbackoption$jobnum -o $logfilepath/matOp_$jobnum.log ../run-larch-usher.sh $larch_usher_options
    done
done

pythonparams=("$inmat" "$iterations" "$numruns" "$logprefix" "+" "${optionnames[@]}")
printf "%s\n" "${pythonparams[@]}" > python_params.txt


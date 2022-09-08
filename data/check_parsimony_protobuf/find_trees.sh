#!/bin/bash
set -eu

OUTDIR=output_trees
FASTA=""
MAX_ALTERNATE_PLACEMENTS=200
DRIFTING_MOVES=4
DRIFTING_TIMES=4
NRUNS=1
VCF=""
REFSEQID=""

print_help()
{
    echo "DESCRIPTION:"
    echo "This script takes a fasta file as input (or a vcf file, and a file \
containing the reference sequence), and uses Usher and matOptimize to \
search for maximally parsimonious trees on those sequences. The output is \
a directory containing many MAT protobufs, each a tree on the same set of \
sequences. The first line of the fasta is expected to contain the name of \
the sequence that will be used as a reference (the sequence of the root \
node of the trees output). \
\
The total number of trees found will be, at maximum, the value passed to '-M' times the value passed to '-D' \
times the value passed to -n."
    echo
    echo "Script requires Usher, matOptimize, and faToVcf."
    echo
    echo "SYNTAX:    find_trees.sh -f INPUT_FASTA [-h|o|M|d|D|v|r]"
    echo
    echo "OPTIONS:"
    echo "-f    Provide an input fasta file (-f or -v REQUIRED)"
    echo "-v    Provide an input vcf file (-f or -v REQUIRED)"
    echo "-r    Provide a fasta file containing a single record: the reference sequence ID and reference sequence (REQUIRED with -v) "
    echo "-n    Specify the number of times to start from scratch rebuilding the tree (default $NRUNS)"
    echo "-o    Specify an output directory for created trees"
    echo "          (default a directory called '$OUTDIR' in the current directory)"
    echo "-M    Specify the maximum number of alternative placements"
    echo "          to be kept when building initial trees. (default $MAX_ALTERNATE_PLACEMENTS)"
    echo "-d    Specify the number of tree moves to apply when drifting (default $DRIFTING_MOVES)"
    echo "-D    Specify the number of times to drift (default $DRIFTING_TIMES)"
    echo "-h    Print this help message and exit"
    echo
}


# getopts expects ':' after options that expect an argument
while getopts "n:f:ho:M:D:d:v:r:" option; do
    case $option in
        n)
            NRUNS=$OPTARG;;
        f)
            FASTA=$OPTARG;;
        h)
            print_help
            exit;;
        o)
            OUTDIR=$OPTARG;;
        M)
            MAX_ALTERNATE_PLACEMENTS=$OPTARG;;
        D)
            DRIFTING_TIMES=$OPTARG;;
        d)
            DRIFTING_MOVES=$OPTARG;;
        v)
            VCF=$OPTARG;;
        r)
            REFSEQID=$OPTARG;;
    esac
done

# TODO: Uncomment this stuff
# [ -e $OUTDIR ] && { echo "$OUTDIR already exists! Exiting."; exit 0; }
mkdir -p $OUTDIR

TMPDIR=$OUTDIR/tmp
mkdir $TMPDIR

if [ -n "${FASTA}" ]
then
    # Sequences are provided in a fasta file:
    REFID=$(head -1 $FASTA)
    REFID=$(echo "${REFID:1}" | xargs)
    VCF=$OUTDIR/out.vcf
    faToVcf $FASTA $VCF -ref=$REFID
else
    # sequences are provided in a vcf file and a reference sequence file:
    [ ! -n "${VCF}" ] && { echo "You must provide an input fasta with '-f' or a vcf with '-v'."; exit 0; }
    [ ! -n "${REFSEQID}" ] && { echo "When using input vcf, you must provide reference sequence ID with '-r'."; exit 0; }
    REFID=$REFSEQID
fi

echo "($REFID)1;" > $TMPDIR/starttree.nh
echo $REFID > $OUTDIR/refid.txt
for ((run=1;run<=NRUNS;run++)); do
    echo Building alternative initial trees: iteration $run / $NRUNS ...
    # place samples in the tree in up to MAX_ALTERNATE_PLACEMENTS different
    # ways
    echo "    Optimizing each resulting tree $DRIFTING_MOVES times ..."
    usher -t $TMPDIR/starttree.nh -v $VCF -o $TMPDIR/mat.pb -d $TMPDIR/ushertree/ -M $MAX_ALTERNATE_PLACEMENTS >& $TMPDIR/usher.log
    for intree in $TMPDIR/ushertree/*.nh; do
        usher -t $intree -v $VCF -o $TMPDIR/mat.pb >> $TMPDIR/usher-optimize.log 2>&1
        for ((optrun=1;optrun<=DRIFTING_TIMES;optrun++)); do
            matOptimize -i $TMPDIR/mat.pb -o $TMPDIR/opt_mat.pb -d $DRIFTING_MOVES >> $TMPDIR/matoptimize.log 2>&1
            mv $TMPDIR/opt_mat.pb $TMPDIR/mat.pb
            cp $TMPDIR/mat.pb $OUTDIR/${run}$(basename $intree)${optrun}.pb
        done
        rm $TMPDIR/mat.pb
        rm -f *intermediate*
    done
    rm -f *intermediate*
    rm -f $TMPDIR/ushertree/*.nh
    rm -f $TMPDIR/ushertree/*.txt
    rm -f $TMPDIR/ushertree/*.tsv
done
rm -r $TMPDIR

#!/bin/bash

# REQUIREMENTS:
# - this script assumes that a working version of the usher executable is accessible in the system's PATH.
# - it also assumes a working python install that has the "Bio" module installed.
# OUTPUTS:
# It creates a directory containing
# - a cleaned fasta file with all gaps written as 'N'
# - a txt file containing a disambiguated consensus sequence for the alignment (computed using the Bio.AlignInfo package)
# - a vcf file corresponding to the cleaned fasta with all IUPAC ambiguity codes overwritten as 'N'
# - a MAT protobuf that is compatible with the vcf and suitable for larch-usher.
if [ $# -lt 2 ]; then
  echo "Usage: [REQUIRED] input fasta file path [REQUIRED] output folder name [OPTIONAL] name of root sequence"
else
  # load the fasta and the directory for outputs
  fasta_file=$1
  output_fp=$2
  # extract the (path-and-extensionless) name of the input fasta file to use as the output filename
  output_fn=$(echo $fasta_file | rev | cut -d "." -f 2 | cut -d "/" -f 1 | rev)
  pb_file="$output_fp/$output_fn.pb"
  vcf_file="$output_fp/$output_fn.vcf"
  mkdir -p $output_fp

  # make sure there are only IUPAC characters in the fasta
  all_ns_fasta_file="$output_fp/$output_fn.fasta"
  sed "s/\-/N/g; s/\?/N/g; s/|/\_/g" $fasta_file > $all_ns_fasta_file

  # can't do this bc the node names are also in the fasta file, and we need to keep those in full
  #sed "s/[^ACGTacgt|]*/N/g; s/|/\_/g" $fasta_file > $all_ns_fasta_file

  if [ $# -lt 3 ]; then 
    # add a consensus sequence to the fasta as a root sequence
    python parse_ambiguous_fasta_for_larch.py $all_ns_fasta_file $output_fp/$output_fn"_reference_sequence"
    root_node="unset_root_seq"
    rs_file=$output_fp/$output_fn"_reference_sequence.txt"
    root_seq=$(head -1 $rs_file)
    sed -i "1s/^/>\ $root_node\n$root_seq\n/" $all_ns_fasta_file
  else
    root_node=$(echo $3 | sed "s/|/\_/; s/\-\?/N/g")
    linenum=$(( $(grep -n "$root_node" $all_ns_fasta_file | cut --delimiter=":" --fields=1)+1 ))
    root_seq=$(awk "NR == $linenum" $all_ns_fasta_file | sed "s/[^ACGTacgt]*/N/g")
    echo $root_seq > "$output_fp/$output_fn_reference_sequence.txt"
  fi

  # output the fasta as a vcf with 'N' instead of more precise ambiguity codes
  faToVcf $all_ns_fasta_file $vcf_file -ambiguousToN
  # now generate a MAT protobuf corresponding to an usher tree constructed on the fasta
  seq1=$(awk 'NR == 1' $all_ns_fasta_file)
  seq2=$(awk 'NR == 3' $all_ns_fasta_file)
  initial_string="(${seq1:2},${seq2:2})1;"
  echo $initial_string > "initial_newick_string_IGNORE.nh"
  USHERCOMMANDS=" -t initial_newick_string_IGNORE.nh -v $vcf_file -o $pb_file -d $output_fp"
  echo $USHERCOMMANDS
  usher $USHERCOMMANDS
fi

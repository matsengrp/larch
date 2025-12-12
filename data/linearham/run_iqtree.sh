#!/bin/sh

# Process each sequence in naive-seqs.fa
# For each: append to igh.fa copy, run iqtree, rename output

naive_fa="naive-seqs.fa"
igh_fa="igh.fa"
tmp_fa="temporary.fa"
n_prefix_len=9

# Generate N prefix string
n_prefix=""
i=0
while [ "$i" -lt "$n_prefix_len" ]; do
    n_prefix="${n_prefix}N"
    i=$((i + 1))
done

seq_name=""
seq_data=""

process_sequence() {
    if [ -z "$seq_name" ]; then
        return
    fi

    # Create temporary file: igh.fa + current sequence (with N prefix)
    cp "$igh_fa" "$tmp_fa"
    printf ">%s\n%s%s\n" "$seq_name" "$n_prefix" "$seq_data" >> "$tmp_fa"

    # Run iqtree
    iqtree -s "$tmp_fa" -o "$seq_name"

    # Rename output treefile to sequence name
    mv "${tmp_fa}.treefile" "${seq_name}.treefile"

    # Clean up iqtree auxiliary files
    rm -f "$tmp_fa" "${tmp_fa}".bionj "${tmp_fa}".ckp.gz "${tmp_fa}".iqtree \
          "${tmp_fa}".log "${tmp_fa}".mldist "${tmp_fa}".model.gz "${tmp_fa}".uniqueseq.phy
}

while IFS= read -r line || [ -n "$line" ]; do
    case "$line" in
        '>'*)
            # Process previous sequence if any
            process_sequence

            # Start new sequence
            seq_name="${line#>}"
            seq_data=""
            ;;
        *)
            # Append to current sequence
            seq_data="${seq_data}${line}"
            ;;
    esac
done < "$naive_fa"

# Process the last sequence
process_sequence

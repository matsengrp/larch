import sys
from Bio import AlignIO
from Bio.Align import AlignInfo

input_fasta_file = sys.argv[1]
output_file_base = sys.argv[2]
output_refseq_file = output_file_base + ".txt"

with open(input_fasta_file, "r") as f:
  algmt = next(AlignIO.parse(f, "fasta"))
  #consensus_ref_seq = AlignInfo.SummaryInfo(algmt).dumb_consensus(ambiguous='N', threshold=0.05/len(algmt))
  consensus_ref_seq = AlignInfo.SummaryInfo(algmt).gap_consensus(ambiguous='N', threshold=0)

print(consensus_ref_seq)
with open(output_refseq_file, "w") as f:
  f.write(str(consensus_ref_seq).replace("N", "A").replace("-", "A"))

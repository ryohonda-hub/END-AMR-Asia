#!/usr/bin/python
#-----------------------
# Check length distribution of a fasta file
#   $ python3 chfa.py <fasta_file>
#-----------------------

import sys
import pandas as pd
from Bio import SeqIO

def main(fasta_file):
    # ======= Check contig length =======
    contig_lengths = [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    if len(contig_lengths) == 0:
        print("No contigs found. Exiting.")
        sys.exit(1)
    contig_lengths.sort(reverse=True)
    # ======= Calculate N50, L50 =======
    total_length = sum(contig_lengths)
    half_length = total_length / 2
    n50 = 0
    l50 = 0
    cumsum = 0
    for idx, length in enumerate(contig_lengths):
        cumsum += length
        if cumsum >= half_length:
            n50 = length
            l50 = idx + 1
            break
    # ======= Calculate average length =======
    mean_length = total_length / len(contig_lengths)
    
    # ======= Create an output dataframe =======
    summary_df = pd.DataFrame({
        'Metric': ['Number of sequences', 'Total length (bp)', 'Mean length (bp)', 'N50 (bp)', 'L50 (count)'],
        'Value': [len(contig_lengths), total_length, round(mean_length, 2), n50, l50]
    })
    # ======= Save the output =======
    print(summary_df.to_string(index=False))
    
if __name__=="__main__":
    args=sys.argv
    main(args[1])

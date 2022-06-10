from sys import argv
from Bio import SeqIO
import pandas as pd

table = pd.read_csv(argv[1], sep='\t', header=None)
table = table[table[0].str.contains('NODE_1')]  # get first contig only
table = table.sort_values([10, 2], ascending=[True, False])
best_ref = table.iloc[0, 1]
with open(argv[2]) as ref, open(argv[3], 'w') as out:
    record_dict = SeqIO.to_dict(SeqIO.parse(argv[2], "fasta"))
    for key in record_dict:
        if best_ref in key:
            SeqIO.write(record_dict[key], out, "fasta")
            break


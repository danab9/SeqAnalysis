import sys
from Bio import SeqIO
import pandas as pd

with open(snakemake.log[0], 'w') as log:
    sys.stderr = sys.stdout = log

    table = pd.read_csv(snakemake.input["table"], sep='\t', header=None)
    table = table[table[0].str.contains('NODE_1')]  # get first contig only
    table = table.sort_values([10, 2], ascending=[True, False])
    best_ref = table.iloc[0, 1]
    with open(snakemake.output["fasta"], 'w') as out:
        record_dict = SeqIO.to_dict(SeqIO.parse(snakemake.input["reference"], "fasta"))
        for key in record_dict:
            if best_ref in key:
                SeqIO.write(record_dict[key], out, "fasta")
                break


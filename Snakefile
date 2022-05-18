configfile : "config/config.yaml"
import pandas as pd

samples = pd.read_csv(config["samples"],index_col="sample", sep ='\t')
IDS=[s for s in list(samples.index)]

include: "rules/bowtie.smk"
include: "rules/samtools.smk"
include: "rules/qc.smk"
include: "rules/msa.smk"

rule all:
    input:
        "msa/alignment.fasta",
        "tree/tree.nwk"
        # expand("fasta/{id}.fa", id=IDS),
        # expand("sam/{id}.sam", id=IDS),
        # expand("qc/fastq/{sample}_fastqc.html",sample=all_fq),
        # expand("qc/trimmed/{sample}_{number}_{paired}_fastqc.html",sample=IDS,number=['1', '2'],paired=['P','UP']),
        # expand("qc/qualimap/{sample}/qualimapReport.html", sample=IDS),



rule rpk:
    input: "stats/{sample}.stats"
    output:
        "stats/{sample}.stats_aug"
    run:
        with open(input[0]) as file:
            with open(output[0],'w') as out:
                for line in file:
                    split_line = line.split("\t")
                    length = int(split_line[1])
                    mapped = int(split_line[2])
                    if length != 0:
                        rpk = (mapped / length) * 1000
                    else:
                        rpk = 0
                    out.write(line.strip() + '\t' + str(rpk) + "\n")





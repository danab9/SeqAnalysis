configfile : "config/config.yaml"
import pandas as pd

samples = pd.read_csv(config["samples"],index_col="sample", sep ='\t')
IDS=[s for s in list(samples.index)]

include: "rules/bowtie.smk"
include: "rules/samtools.smk"
include: "rules/qc.smk"
include: "rules/msa.smk"
include: "rules/denovo.smk"
include: "rules/artificialref.smk"
include: "rules/screen.smk"
include: "rules/decontamination.smk"

rule all:
    input:
        #"msa/alignment.fasta",
        #"tree/tree.nwk",
        #"variability/variability.txt",
        #"variability/variability.png",
        #"tree/tree.png",
        #"reference/artificial_reference_ERR4082860.fa",
        expand("reference/artificial_reference_{sample}.fa", sample=IDS),
        expand("best_references/{sample}.fasta", sample=IDS),
        "qc/multiqc_report.html" #als add




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





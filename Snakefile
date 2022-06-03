configfile : "config/config.yaml"
import pandas as pd

samples = pd.read_csv(config["samples"],index_col="sample", sep ='\t')
IDS=[s for s in list(samples.index)]

include: "rules/bowtie.smk"
include: "rules/samtools.smk"
include: "rules/qc.smk"
include: "rules/msa.smk"
include: "rules/denovo.smk"

rule all:
    input:
        #"msa/alignment.fasta",
        #"tree/tree.nwk",
        #"variability/variability.txt",
        #"variability/variability.png",
        #"tree/tree.png",
        expand("denovo_assembly/{sample}/contigs.fasta", sample=IDS)


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





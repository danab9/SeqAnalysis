configfile : "config/config.yaml"
import pandas as pd


# samples = pd.read_csv(config["samples"], index_col = "sample", sep ='\t')
samples = pd.read_csv("samples.tsv",index_col=0, sep =',')
IDS=[x+"_tiny" for x in list(samples.index)]
r1 = lambda wildcards:samples.at[wildcards.sample, 'fq1']
r2 = lambda wildcards:samples.at[wildcards.sample, 'fq2']  # for later
print(r1)

fastq_dir = "fastq"
ref_prefix = config['ref']

include: "rules/bowtie.smk"
include: "rules/samtools.smk"

rule all:
    input: expand("results/{id}.rpk", id=IDS)

rule rpk:
    input: "stats/{dana}.stats"
    output:
        "results/{dana}.rpk"
    run:
        with open(input[0]) as file:
            with open(output[0],'w') as out:
                for line in file:
                    split_line = line.split("\t")
                    length = int(split_line[1])
                    mapped = int(split_line[2])
                    if length != 0:
                        rpk = (mapped / length) / 1000
                    else:
                        rpk = 0
                    out.write(line + '\t' + str(rpk))


    #     input: expand("sam_tiny/{id}.sam", id=IDS)


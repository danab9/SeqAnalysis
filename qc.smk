configfile : "config/config.yaml"

import pandas as pd
samples = pd.read_csv(config["samples"],index_col="sample", sep ='\t')
IDS=[s for s in list(samples.index)]
r1 = lambda wildcards:samples.at[wildcards.sample, 'fq1']
r2 = lambda wildcards:samples.at[wildcards.sample, 'fq2']
all_fq = [ID + "_tiny_1" for ID in IDS] + [ID + "_tiny_2" for ID in IDS]  # todo!


rule temp:
    input:
        qc=expand("fastqc/{sample}_fastqc.html", sample=all_fq),
        trimmed_r1_p= expand("trimmed/{myrthe}_1_P.fastq.gz", myrthe=IDS)


rule fastqc:
    input:
        fq="fastq/{sample}.fastq.gz"
    output:
        html="fastqc/{sample}_fastqc.html",
        zip="fastqc/{sample}_fastqc.zip"
    conda:
        config["env"]
    log:
        "logs/fastqc/{sample}.log"
    shell:
        "fastqc {input.fq} -o fastqc &> {log}"


rule trimmomatic:
    input:
        r1 = lambda wildcards: samples.at[wildcards.sample, 'fq1'],
        r2 = lambda wildcards: samples.at[wildcards.sample, 'fq2']
    output:
        r1_p= "trimmed/{sample}_1_P.fastq.gz", r1_u = "trimmed/{sample}_1_UP.fastq.gz",
        r2_p= "trimmed/{sample}_2_P.fastq.gz", r2_u = "trimmed/{sample}_2_UP.fastq.gz"
    # TODO: add params (through config + params)
    conda:
        config["env"]
    log:
        "logs/trimmomatic/{sample}.log"
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1_p} {output.r1_u} {output.r2_p} {output.r2_u} TRAILING:20 &> {log}"


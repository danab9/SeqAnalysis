configfile : "config/config.yaml"

import pandas as pd
samples = pd.read_csv(config["samples"],index_col="sample", sep ='\t')
IDS=[s for s in list(samples.index)]
r1 = lambda wildcards:samples.at[wildcards.sample, 'fq1']
r2 = lambda wildcards:samples.at[wildcards.sample, 'fq2']
all_fq = [ID + "_1" for ID in IDS] + [ID + "_2" for ID in IDS]  # todo!
env_qc = "../envs/multiqc.yaml"

rule temp:
    input:
        qc=expand("qc/fastq/{sample}_fastqc.html", sample=all_fq),
        trimmed_r1_p= expand("trimmed/{myrthe}_1_P.fastq.gz", myrthe=IDS),
        trimmed_qc=expand("qc/trimmed/{sample}_{number}_{paired}_fastqc.html", sample=IDS, number=['1','2'], paired=['P','UP'])


rule fastqc:
    input:
        fq="{folder}/{sample}.fastq.gz"
    output:
        html="qc/{folder}/{sample}_fastqc.html",
        zip="qc/{folder}/{sample}_fastqc.zip"
    conda:
        config["env"]
    log:
        "logs/qc/{folder}/{sample}.log"
    shell:
        "fastqc {input.fq} -o qc/{wildcards.folder} &> {log}"


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


rule qualimap:
    input:
        "bam_sorted/{sample}_sorted.bam"
    output:
        dir=directory("qc/qualimap/{sample}"),
        file=("qc/qualimap/{sample}/qualimapReport.html"),
    conda:
        config["env"]
    log:
        "logs/qualimap/{sample}.log"
    shell:
        "qualimap bamqc -bam {input} -outdir {output.dir} &> {log}"

rule multiqc:
    input:
        fqc=expand("qc/fastq/{sample}_fastqc.html",sample=all_fq) if config["skip_trimming"] in['False',''] else [],
        tqc=expand("qc/trimmed/{sample}_{number}_{paired}_fastqc.html",sample=IDS,number=['1', '2'],paired=['P', 'UP'])
            if config['skip_fastQC'] in ['False',''] else [],
        qualimap=expand("qc/qualimap/{sample}/qualimapReport.html",sample=IDS) if config['skip_qualimap'] in
                                                                                  ['False',''] else []
    output:
        "qc/multiqc_report.html"
    conda:
        config["multiqc_env"]
    log:
        "logs/multiqc/multiqc.log"
    params:
        config["multiqcparam"]  # for example: -f parameter to ensure existing multiqc report is override.
    shell:
        "multiqc qc {params} -o qc &> {log}"

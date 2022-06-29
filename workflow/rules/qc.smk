rule fastqc_untrimmed:
    input:
        fq = lambda wildcards: samples.at[wildcards.sample, 'fq'+ wildcards.number],
    output:
        html="results/qc/fastq/{sample}_{number}_fastqc.html",
        zip="results/qc/fastq/{sample}_{number}_fastqc.zip"
    conda:
        "../envs/env.yaml"
    log:
        "results/logs/qc/untrimmed/{sample}_{number}.log"
    shell:
        "fastqc {input.fq} -o results/qc/fastq &> {log}"

rule fastqc_trimmed:
    input:
        "results/fastq/trimmed/{sample}_{number}_{paired}.fastq.gz"
    output:
        html="results/qc/trimmed/{sample}_{number}_{paired}_fastqc.html",
        zip="results/qc/trimmed/{sample}_{number}_{paired}_fastqc.zip"
    conda:
        "../envs/env.yaml"
    log:
        "results/logs/qc/trimmed/{sample}_{number}_{paired}.log"
    shell:
        "fastqc {input} -o results/qc/trimmed &> {log}"


rule trimmomatic:
    input:
        r1 = lambda wildcards: samples.at[wildcards.sample, 'fq1'],
        r2 = lambda wildcards: samples.at[wildcards.sample, 'fq2']
    output:
        r1_p= "results/fastq/trimmed/{sample}_1_P.fastq.gz", r1_u = "results/fastq/trimmed/{sample}_1_UP.fastq.gz",
        r2_p= "results/fastq/trimmed/{sample}_2_P.fastq.gz", r2_u = "results/fastq/trimmed/{sample}_2_UP.fastq.gz"
    params:
        trailing = config["trimmomatic"]['trailing'],
        illuminaclip = ':'.join(config["trimmomatic"]["illuminaclip"].values())
    conda:
        "../envs/env.yaml"
    log:
        "results/logs/trimmomatic/{sample}.log"
    threads: 4
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1_p} {output.r1_u} {output.r2_p} {output.r2_u} TRAILING:{params.trailing} ILLUMINACLIP:{params.illuminaclip} -threads {threads} &> {log}"


rule qualimap:
    input:
        "results/bam/sorted/{sample}.bam"
    output:
        dir=directory("results/qc/qualimap/{sample}"),
        file="results/qc/qualimap/{sample}/qualimapReport.html"
    conda:
        "../envs/env.yaml"
    log:
        "results/logs/qualimap/{sample}.log"
    threads: 4
    shell:
        "qualimap bamqc -bam {input} -outdir {output.dir} -nt {threads} &> {log}"

rule multiqc:
    input:
        fqc=expand("results/qc/fastq/{sample}_{number}_fastqc.html",sample=IDS,number=['1', '2']) if config["skip_trimming"] in['False',''] else [],
        tqc=expand("results/qc/trimmed/{sample}_{number}_{paired}_fastqc.html",sample=IDS,number=['1', '2'],paired=['P', 'UP'])
            if config['skip_fastQC'] in ['False',''] else [],
        qualimap=expand("results/qc/qualimap/{sample}/qualimapReport.html",sample=IDS) if config['skip_qualimap'] in
                                                                                  ['False',''] else []
    output:
        "results/qc/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    log:
        "results/logs/multiqc/multiqc.log"
    threads: 1
    params:
        config["multiqcparam"]
    shell:
        "multiqc results/qc {params} -o results/qc &> {log}" # fixed 1 mistake that we did not submit previously, which probably caused the output to the shell.

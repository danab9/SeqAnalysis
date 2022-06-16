r1 = lambda wildcards:samples.at[wildcards.sample, 'fq1']
r2 = lambda wildcards:samples.at[wildcards.sample, 'fq2']


rule temp:  # what's that?
    input:
        qc=expand("../report/qc/fastq/{sample}_fastqc.html", sample=all_fq),
        trimmed_r1_p= expand("../results/trimmed/{myrthe}_1_P.fastq.gz", myrthe=IDS),
        trimmed_qc=expand("../report/qc/trimmed/{sample}_{number}_{paired}_fastqc.html", sample=IDS, number=['1','2'], paired=['P','UP'])


rule fastqc:
    input:
        fq="{folder}/{sample}.fastq.gz"
    output:
        html="../report/qc/{folder}/{sample}_fastqc.html",
        zip="../report/qc/{folder}/{sample}_fastqc.zip"
    conda:
        "../envs/env.yaml"
    log:
        "../results/logs/qc/{folder}/{sample}.log"
    shell:
        "fastqc {input.fq} -o qc/{wildcards.folder} &> {log}"


rule trimmomatic:
    input:
        r1 = lambda wildcards: samples.at[wildcards.sample, 'fq1'],
        r2 = lambda wildcards: samples.at[wildcards.sample, 'fq2']
    output:
        r1_p= "../results/fastq/trimmed/{sample}_1_P.fastq.gz", r1_u = "../results/fastq/trimmed/{sample}_1_UP.fastq.gz",
        r2_p= "../results/fastq/trimmed/{sample}_2_P.fastq.gz", r2_u = "../results/fastq/trimmed/{sample}_2_UP.fastq.gz"
    params:
        trailing = config["trimmomatic"]['trailing'],
        illuminaclip = ':'.join(config["trimmomatic"]["illuminaclip"].values())
    conda:
        "../envs/env.yaml"
    log:
        "../results/logs/trimmomatic/{sample}.log"
    threads: 4
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1_p} {output.r1_u} {output.r2_p} {output.r2_u} TRAILING:{params.trailing} ILLUMINACLIP:{params.illuminaclip} -threads {threads} &> {log}"


rule qualimap:
    input:
        "../results/bam/sorted/{sample}.bam"
    output:
        dir=directory("../report/qc/qualimap/{sample}"),
        file="../report/qc/qualimap/{sample}/qualimapReport.html"
    conda:
        "../envs/env.yaml"
    log:
        "../results/logs/qualimap/{sample}.log"
    threads: 4
    shell:
        "qualimap bamqc -bam {input} -outdir {output.dir} -nt {threads} &> {log}"

rule multiqc:
    input:
        fqc=expand("../report/qc/fastq/{sample}_fastqc.html",sample=all_fq) if config["skip_trimming"] in['False',''] else [],
        tqc=expand("../report/qc/trimmed/{sample}_{number}_{paired}_fastqc.html",sample=IDS,number=['1', '2'],paired=['P', 'UP'])
            if config['skip_fastQC'] in ['False',''] else [],
        qualimap=expand("../report/qc/qualimap/{sample}/qualimapReport.html",sample=IDS) if config['skip_qualimap'] in
                                                                                  ['False',''] else []
    output:
        "../report/qc/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    log:
        "../results/logs/multiqc/multiqc.log"
    threads: 1
    params:
        config["multiqcparam"]  # for example: -f parameter to ensure existing multiqc report is override.
    shell:
        "multiqc qc {params} -o qc &> {log}"

rule kraken:
    input:
        r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1']) if config[
            "skip_trimming"]=='True' else "../results/fastq/trimmed/{sample}_1_P.fastq.gz",
        r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2']) if config[
            "skip_trimming"]=='True' else "../results/fastq/trimmed/{sample}_2_P.fastq.gz"
    output:
        clasified_reads_1 = "../results/fastq/kraken/{sample}_1.fq",
        clasified_reads_2 = "../results/fastq/kraken/{sample}_2.fq",
        report = "report/kraken/{sample}.txt"
    log:
        "../results/logs/kraken/{sample}.log"
    threads: 10
    params:
        kraken_db = config["kraken_db"], #Cannot use as input, because; The flag 'directory' used in rule screening is only valid for outputs, not inputs.
    conda:
        "../envs/decontamination.yaml"
    shell:
        "kraken2 -db {params.kraken_db} --paired --classified-out ../results/fastq/kraken/{wildcards.sample}#.fq {input.r1} {input.r2} --report {output.report} --threads {threads} &> {log}"

rule screen: #creates a multiqc report of the screened reads in the folder qc/screened
    input:
        screen = expand("report/kraken/{sample}.txt",sample=IDS),
    output:
        "report/kraken/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    log:
        "../results/logs/multiqc/multiqc_kraken.log"
    threads: 1
    params:
        config["multiqcparam"]  # for example: -f parameter to ensure existing multiqc report is override.
    shell:
        "multiqc report/kraken {params} -o report/kraken &> {log}"



# rule test:
#     input:
#         config["references"]
#     output:
#         "../results/test2.txt"
#     conda:
#         "../envs/blast.yaml" #because we should not use variables for environments
#     log:
#         "../results/logs/test.log"
#     shell:
#         "touch {output}"

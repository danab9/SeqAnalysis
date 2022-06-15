rule kraken:
    input:
        r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1']) if config[
            "skip_trimming"] else "trimmed/{sample}_1_P.fastq.gz",
        r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2']) if config[
            "skip_trimming"] else "trimmed/{sample}_2_P.fastq.gz"
    output:
        clasified_reads_1 = "screened/{sample}_1_P.fq",
        clasified_reads_2= "screened/{sample}_2_P.fq",
        report = "screened/{sample}.txt"
    log:
        "logs/kraken/{sample}.log"
    threads: 10
    params:
        kraken_db = config["kraken_db"], #Cannot use as input, because; The flag 'directory' used in rule screening is only valid for outputs, not inputs.
    conda:
        "../envs/decontamination.yaml"
    shell:
        "kraken2 -db {params.kraken_db} --paired --classified-out screened/{wildcards.sample}#_P.fq {input.r1} {input.r2} --report {output.report} --threads {threads} &> {log}"

rule screen: #creates a multiqc report of the screened reads in the folder qc/screened
    input:
        screen=expand("screened/{sample}.txt",sample=IDS),
    output:
        "qc/screened/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/multiqc/multiqc_screen.log"
    threads: 1
    params:
        config["multiqcparam"]  # for example: -f parameter to ensure existing multiqc report is override.
    shell:
        "multiqc screened {params} -o qc/screened &> {log}"

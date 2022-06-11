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

rule screen:
    input:
        screen=expand("screened/{sample}.txt",sample=all_fq) if config["decontamination"] in['False',''] else [],
    output:
        "qc/multiqc_krakenscreen_report.html"
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/multiqc/multiqc.log"
    threads: 1
    params:
        config["multiqcparam"]  # for example: -f parameter to ensure existing multiqc report is override.
    shell:
        "multiqc qc {params} -o qc &> {log}"
# https://telatin.github.io/microbiome-bioinformatics/MultiQC/
# https://github.com/DerrickWood/kraken2/wiki/Manual#classification use --report

# rule contamination_fasta:
#     #create a fasta from the output report by kraken2, for each sample.
#     # In the assignmetn is not stated that we have to do this, but that the user should do this based on the MULTIQC
#

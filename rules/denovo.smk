rule spades:  # denovo assembly
    input:
        r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1']) if config[
        "skip_trimming"] and not config["decontamination"] == "True" else "screened_fastq/{sample}_screened_1.fq"
        if config["decontamination"] == "True" else "trimmed/{sample}_1_P.fastq.gz",
        r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2']) if config[
        "skip_trimming"] and not config["decontamination"] == "True" else "screened_fastq/{sample}_screened_2.fq'"
        if config["decontamination"] == "True" else "trimmed/{sample}_2_P.fastq.gz"
    output:
        file="denovo_assembly/{sample}/contigs.fasta"
    conda:
        "../envs/spades.yaml"
    log:
        "logs/spades/{sample}.log"
    threads: 4
    shell:
        "spades.py -1 {input.r1} -2 {input.r2} -o denovo_assembly/{wildcards.sample} --threads {threads} &>{log}"



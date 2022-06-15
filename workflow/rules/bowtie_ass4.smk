rule bowtie2_build:
    input:
        "reference/NC_045512.2.fa"
    output:
        multiext(
            "reference/NC_045512.2",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )
    log:
        "logs/bowtie2_build/build_{sample}.log"
    threads: 4
    conda:
        "../envs/env.yaml"
    shell:
        "bowtie2-build {input} --threads {threads} &> {log}"

rule mapreads:
    input:
        indexed_ref = multiext(
                "reference/artificial_reference_{sample}",
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2"),
        r1 = (lambda wildcards:samples.at[wildcards.sample, 'fq1']) if config[
            "skip_trimming"] else "trimmed/{sample}_1_P.fastq.gz",
        r2 = (lambda wildcards:samples.at[wildcards.sample, 'fq2']) if config[
            "skip_trimming"] else "trimmed/{sample}_2_P.fastq.gz"
    output:
        "sam/{sample}.sam"
    params:
        local = config["bowtie2"]["local"],
        ma = config["bowtie2"]["ma"]
    log:
        "logs/bowtie2/{sample}_aligment.log"
    threads: 6
    conda:
        "../envs/env.yaml"
    shell:
        "bowtie2 -x reference/artificial_reference_{wildcards.sample} -1 {input.r1} -2 {input.r2} -S {output} --threads {threads} &> {log}"

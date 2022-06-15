rule bowtie2_build:
    input:
        "../results/references/artificial/consensus/{sample}.fa"
    output:
        multiext(
            "../results/references/artificial/consensus/{sample}",  # .fa?
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )
    log:
        "../results/logs/bowtie2_build/build_{sample}.log"
    threads: 4
    conda:
        "../envs/env.yaml"
    shell:
        "bowtie2-build {input} ../results/reference/artificial/consensus/{wildcards.sample} --threads {threads} &> {log}"

rule mapreads:
    input:
        indexed_ref = multiext(
                "../results/references/artificial/consensus/{sample}",
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2"),
        r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1']) if config[
            "skip_trimming"] and not config["decontamination"]=="True" else "../results/fastq/decontaminated/{sample}_1.fq"
            if config["decontamination"]=="True"  else "../results/fastq/trimmed/{sample}_1_P.fastq.gz",
        r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2']) if config[
            "skip_trimming"] and not config["decontamination"]=="True" else "../results/fastq/decontaminated/{sample}_2.fq"
            if config["decontamination"]=="True" else "../results/fastq/trimmed/{sample}_2_P.fastq.gz" #"decontaminated/{sample}_2_P.fastq.gz" if config[decontamination] elseif .. else ..
    output:
        "../results/sam/{sample}.sam"
    params:
        local = config["bowtie2"]["local"],
        ma = config["bowtie2"]["ma"]
    log:
        "../results/logs/bowtie2/{sample}_aligment.log"
    threads: 6
    conda:
        "../envs/env.yaml"
    shell:
        "bowtie2 -x ../results/references/artificial/consensus/{sample} -1 {input.r1} -2 {input.r2} -S {output} --ma {params.ma} {params.local} --threads {threads} &> {log}"
# bowtie

rule bowtie_map_contaminations:
    input:
        indexed_ref= multiext(
            "reference/to_screen.fa",  # TODO: which reference is removed?
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2"),
        r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1']) if config[
                "skip_trimming"] else "trimmed/{sample}_1_P.fastq.gz",
        r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2']) if config[
                "skip_trimming"] else "trimmed/{sample}_2_P.fastq.gz"  #"decontaminated/{sample}_2_P.fastq.gz" if config[decontamination] elseif .. else ..
    output:
        "sam_contaminations/{sample}.sam"
    params:
        local=config["bowtie2"]["local"],
        ma=config["bowtie2"]["ma"]
    log:
        "logs/bowtie2/{sample}_screen_aligment.log"
    threads: 6
    conda:
        "../envs/env.yaml"
    shell:
        "bowtie2 -x reference/to_screen.fa -1 {input.r1} -2 {input.r2} -S {output} --threads {threads} &> {log}"

# samtools
rule keep_unmapped:   # TODO: add to rule all, see how to use bowtie mapping rule differently each time
    input:
        "sam_contaminations/{sample}.sam"
    output:
        "bam/{sample}_unmapped.bam"
    conda:
        "../envs/env.yaml"
    threads: 4
    log:
        "logs/samtools/{sample}_view_unmapped.log"
    shell:
        "samtools view -b -f 4 {input} --threads {threads} > {output} 2> {log}"

rule sam_to_fastq:
    input:
        "bam/{sample}_unmapped.bam"
    output:
        fq1="fastq/{sample}_screened_1.fq", fq2="fastq/{sample}_screened_2.fq"
    conda:
        "../envs/decontamination.smk"
    log:
        "logs/bamtofq/{sample}.log"
    shell:
        "bedtoold bamtofastq -i {input} -fq {output.fq1} -fq2 {output.fq2} 2> {log}"
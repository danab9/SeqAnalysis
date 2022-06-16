# bowtie
to_screen_ref=config["contamination_reference"]

rule bowtie2_build_contamination:
    input:
        config["contamination_reference"]
    output:
        multiext(
            config["contamination_reference"],
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )
    log:
        "../results/logs/bowtie2_build/build_contamination.log"
    threads: 4
    conda:
        "../envs/env.yaml"
    shell:
        "bowtie2-build {input} --threads {threads} &> {log}"

rule bowtie_map_contaminations:
    input:
        ref=to_screen_ref,
        indexed_ref= multiext(
            to_screen_ref,
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2"),
        r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1']) if config[
                "skip_trimming"]=='True' else "../results/fastq/trimmed/{sample}_1_P.fastq.gz",
        r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2']) if config[
                "skip_trimming"]=='True' else "../results/fastq/trimmed/{sample}_2_P.fastq.gz"  #"decontaminated/{sample}_2_P.fastq.gz" if config[decontamination] elseif .. else ..
    output:
        "../results/sam/contaminations/{sample}.sam"
    params:
        local=config["bowtie2"]["local"],
        ma=config["bowtie2"]["ma"]
    log:
        "../results/logs/bowtie2/contamination_alignment/{sample}.log"
    threads: 6
    conda:
        "../envs/env.yaml"
    shell:
        "bowtie2 -x {input.ref} -1 {input.r1} -2 {input.r2} -S {output} --threads {threads} &> {log}"

# samtools
rule keep_unmapped:   # TODO: add to rule all, see how to use bowtie mapping rule differently each time
    input:
        "../results/sam/contaminations/{sample}.sam"
    output:
        "../results/bam/contaminations/{sample}_unmapped.bam"
    conda:
        "../envs/env.yaml"
    threads: 4
    log:
        "../results/logs/samtools/contaminations/{sample}_unmapped.log"
    shell:
        "samtools view -b -f 4 {input} --threads {threads} > {output} 2> {log}"

rule sam_to_fastq:
    input:
         "../results/bam/contaminations/{sample}_unmapped.bam"
    output:
        fq1="../results/fastq/decontaminated/{sample}_1.fq", fq2="../results/fastq/decontaminated/{sample}_2.fq"
    conda:
        "../envs/decontamination.smk"
    log:
        "../results/logs/bamtofq/{sample}_decontaminated.log"
    shell:
        "bedtoold bamtofastq -i {input} -fq {output.fq1} -fq2 {output.fq2} 2> {log}"
rule bowtie2_build_contamination:
    input:
        config["contamination_reference"]
    output:
        multiext(
            "results/references/contamination/contamination_reference",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )
    log:
        "results/logs/bowtie2_build/build_contamination.log"
    threads: 4
    conda:
        "../envs/env.yaml"
    shell:
        """
        mkdir -p results/references/contamination && cp {input} results/references/contamination/contamination_reference.fa
        bowtie2-build {input} results/references/contamination/contamination_reference --threads {threads} &> {log}
        """
    # Error No output file specified!, to fix this I added {input} for a 2nd time

rule bowtie_map_contaminations:
    input:
        ref="results/references/contamination/contamination_reference.fa",
        #@Sandro: It is true that we don't use the config defined directory here, but this dir is used in the build rule above,
        # this rule copies the user defined file name to this new location, so that bowtie build can use it.
        indexed_ref= multiext(
            "results/references/contamination/contamination_reference",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2"),
        r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1'])
                if config["skip_trimming"]=='True' else "results/fastq/trimmed/{sample}_1_P.fastq.gz",
        r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2'])
                if config["skip_trimming"]=='True' else "results/fastq/trimmed/{sample}_2_P.fastq.gz"
    output:
        "results/sam_contaminations/{sample}.sam"
    log:
        "results/logs/bowtie2/contamination_alignment/{sample}.log"
    threads: 6
    conda:
        "../envs/env.yaml"
    shell:
        "bowtie2 -x results/references/contamination/contamination_reference -1 {input.r1} -2 {input.r2} -S {output} --threads {threads} &> {log}"

rule keep_unmapped:   # TODO: add to rule all, see how to use bowtie mapping rule differently each time
    input:
        "results/sam_contaminations/{sample}.sam"
    output:
        "results/bam_decontaminated/{sample}.bam"
    conda:
        "../envs/env.yaml"
    threads: 4
    log:
        "results/logs/samtools/contaminations/{sample}_bam_unmapped.log"
    shell:
        "samtools view -b -f 4 {input} --threads {threads} > {output} 2> {log}"

rule sam_to_fastq:
    input:
         "results/bam_decontaminated/{sample}.bam"
    output:
        fq1="results/fastq/decontaminated/{sample}_1.fq", fq2="results/fastq/decontaminated/{sample}_2.fq"
    conda:
        "../envs/env.yaml"
    log:
        "results/logs/bamtofq/{sample}_decontaminated.log"
    shell:
        "bedtools bamtofastq -i {input} -fq {output.fq1} -fq2 {output.fq2} 2> {log}"
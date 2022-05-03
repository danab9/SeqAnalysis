import pandas as pd
sam_dir = "sam_tiny"
configfile : "config/config.yaml"
samples = pd.read_csv("samples.tsv",index_col=0, sep =',')
IDS=[x+"_tiny" for x in list(samples.index)]
# r1 = lambda wildcards:samples.at[wildcards.sample, 'fq1']
fastq_dir = "fastq"
ref_prefix = config['ref']

rule samtobam:
    input:
        sam_dir+"/{dana}.sam"
    output:
        "bam/{dana}.bam"
    conda:
       "envs/env.yaml"
    threads: 4
    shell:
      "env samtools view -b {input} > {output}"

# rule samtobam:
#     input:
#         sam_dir+"/{dana}.sam"
#     output:
#         "bam/{dana}.bam"
#     shell:
#         "samtools view -b {input} > {output}"

rule sort:
    input:
        "bam/{dana}.bam"
    output:
        "bam_sorted/{dana}.bam.sorted"
    threads: 4
    shell:
        "samtools sort {input} -o {output}"

rule index:
    input:
        "bam_sorted/{dana}.bam.sorted"
    output:
        "bam_sorted/{dana}.bam.sorted.bai"
    threads: 4
    shell:
        "samtools index -b {input}"

rule mapstats:
    input:
        sorted="bam_sorted/{dana}.bam.sorted",
        index="bam_sorted/{dana}.bam.sorted.bai"
    output:
        "stats/{dana}.stats"
    threads: 4
    shell:
        "samtools idxstats {input.sorted} > {output}"

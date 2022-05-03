import pandas as pd
configfile : "config/config.yaml"

samples = pd.read_csv("samples.tsv",index_col=0, sep =',')
IDS=[x+"_tiny" for x in list(samples.index)]
# r1 = lambda wildcards:samples.at[wildcards.sample, 'fq1']
fastq_dir = "fastq"
ref_prefix = config['ref']
sam_dir = "sam_tiny"


rule bowtie2_build:
    input:
        ref=ref_prefix+".fa"
    output:
        multiext(
            ref_prefix,
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        "logs/bowtie2_build/build.log",
    params:
        extra="",  # optional parameters
    threads: 4
    wrapper:
        "v1.3.2/bio/bowtie2/build"

rule mapreads:
    input:
        indexed_ref = multiext(
                ref_prefix,
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2"),
        read1 = fastq_dir + '/{dana}_1.fastq.gz',
        read2 = fastq_dir + '/{dana}_2.fastq.gz'
    output:
        sam_dir+"/{dana}.sam"
    threads: 4
    shell:
        "bowtie2 -x {ref_prefix} -1 {input.read1} -2 {input.read2} -S {output}"
import pandas as pd
configfile : "config/config.yaml"

samples = pd.read_csv(config["samples"],index_col="sample", sep ='\t')
ref_prefix = config['ref']

rule bowtie2_build:
    input:
        ref_prefix+".fa"
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
        "logs/bowtie2_build/build.log"
    threads: 4
    conda:
        config["env"]
    shell:
        "bowtie2-build {input} {ref_prefix} --threads {threads} &> {log}"

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
        r1 = lambda wildcards:samples.at[wildcards.sample, 'fq1'],
        r2 = lambda wildcards:samples.at[wildcards.sample, 'fq2']
    output:
        "sam/{sample}.sam"
    log:
        "logs/bowtie2/{sample}_aligment.log"
    threads: 4
    conda:
        config["env"]
    shell:
        "bowtie2 -x {ref_prefix} -1 {input.r1} -2 {input.r2} -S {output} --threads {threads} &> {log}"


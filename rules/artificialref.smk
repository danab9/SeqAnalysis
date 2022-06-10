import pandas as pd
configfile : "config/config.yaml"
blast_env = "../envs/blast.yaml"

samples = pd.read_csv(config["samples"],index_col="sample", sep ='\t')

rule mapcontigs:
    input:
        references = "reference.fa", #fasta file with multiple user provided references
        contigs = "denovo_assembly/{sample}/contigs.fasta"
    output:
        best_reference = "best_reference_{sample}.fa"
    log:
        "logs/mapcontigs.log"
    threads: 4
    conda:
        blast_env
    shell:
        #BLAST


import pandas as pd
configfile : "config/config.yaml"
samples = pd.read_csv(config["samples"],index_col="sample", sep ='\t')
ref_prefix = config['ref']
env = "../envs/artificialref.yaml"

rule readblast:
    input:
        references = "reference.fa", #fasta file with multiple user provided references
        contigs = "denovo_assembly/{sample}/contigs.fasta"
    output:
        best_reference = "best_reference_{sample}.fa"
    log:
        "logs/mapcontigs.log"
    threads: 4
    conda:
        env
    shell:
        #

rule artificialreference:
    input: # contigs for the sample, and best reference for the sample
        best_reference = "best_reference_{sample}.fa", # reference that was most similar to our sample.
        contigs = "denovo_assembly/{sample}/contigs.fasta"
    output:
        artificial_reference = "artificial_reference_{sample}.sam"
    log:
        "logs/artificialreference.log"
    threads: 4
    conda:
        env
    shell:
        "minimap2 -a {input.best_reference} {input.contigs} > {output.artificial_reference} 2> {log}"
        # minimap https://github.com/lh3/minimap2

rule artificialrefconcensus:
    input:
        best_reference = "best_reference_{sample}.fa",
        sam = "artificial_reference_{sample}.sam"
    output:
        consensus = "artificial_reference_{sample}.fa"
    log:
        "logs/bcsf.log"
    threads: 4
    conda:
        env
    shell:
        "bcftools mpileup -B -Ou -f {input.best_reference} {input.sam} | bcftools call -mv -M -Oz -o calls.vcf.gz 2> {log}"
        "bcftools index calls.vcf.gz 2> {log}"
        "cat {input.best_reference} | bcftools consensus calls.vcf.gz > {output.consensus} 2> {log}"
        # TODO; check if it works with 2logs and multiple commands.
import os

ref_name = "results/references/" + os.path.basename(config["references"]).split(".")[0]
rule makeblastdb:
    input:
        ref_path = config["references"]
    output:
        multiext(
            ref_name,
            ".nhr",
            ".nin",
            ".nog",
            ".nsd",
            ".nsi",
            ".nsq"
        )
    conda:
        "../envs/blast.yaml"
    log:
        "results/logs/blastdb.log"
    shell:
        """
        makeblastdb -in {input.ref_path} -out {ref_name} -dbtype nucl -parse_seqids -logfile {log}
        """

rule mapcontigs:
    input:
        contigs = "results/denovo_assembly/{sample}/contigs.fasta",
        indexed_ref = multiext(
            ref_name,
            ".nhr",
            ".nin",
            ".nog",
            ".nsd",
            ".nsi",
            ".nsq"
        )
    output:
        "results/blast/{sample}.tsv"
    log:
        "results/logs/blast/{sample}.log"
    threads: 4
    conda:
        "../envs/blast.yaml"
    shell:
        "blastn -query {input.contigs} -db {ref_name} -outfmt 6 -out {output} -num_threads {threads} 2>{log}"


rule best_reference:
    input:
        table = "results/blast/{sample}.tsv",
        reference = config["references"]
    output:
        fasta = "results/references/best_references/{sample}.fasta"
    conda:
        "../envs/artificialref.yaml"
    log:
        "results/logs/best_references/{sample}.log"
    script:
        "../scripts/best_reference.py"


rule artificial_reference:
    input:
        best_reference = "results/references/best_references/{sample}.fasta",
        contigs = "results/denovo_assembly/{sample}/contigs.fasta"
    output:
        artificial_reference = "results/references/artificial/bam/{sample}.bam"
    log:
        "results/logs/artificialreference/{sample}.log"
    threads: 4
    conda:
        "../envs/artificialref.yaml"
    shell:
        "minimap2 -t {threads} -a {input.best_reference} {input.contigs} 2> {log} | samtools view -b - > {output.artificial_reference} 2> {log}"

rule sort_artificial_reference:
    input:
        "results/references/artificial/bam/{sample}.bam"
    output:
        "results/references/artificial/bam_sorted/{sample}.bam"
    log:
        "results/logs/samtools/artificial/{sample}_sort.log"
    threads: 4
    conda:
        "../envs/env.yaml"
    shell:
        "samtools sort {input} -o {output} --threads {threads} &> {log}"


rule artificialrefconsensus:
    input:
        best_reference = "results/references/best_references/{sample}.fasta",
        bam_sorted = "results/references/artificial/bam_sorted/{sample}.bam"
    output:
        consensus = "results/references/artificial/{sample}.fa"
    log:
        "results/logs/bcsf/{sample}.log"
    threads: 1
    conda:
        "../envs/artificialref.yaml"
    shell:
        """
        bcftools mpileup -B -Ou -f {input.best_reference} {input.bam_sorted} | bcftools call -mv -M -Oz -o results/{wildcards.sample}_calls.vcf.gz 2> {log}
        bcftools index results/{wildcards.sample}_calls.vcf.gz -f 2> {log}
        cat {input.best_reference} | bcftools consensus results/{wildcards.sample}_calls.vcf.gz > {output.consensus} 2> {log}
        """


rule makeblastdb:
    input:
        "reference/reference.fa"
    output:
        multiext(
            'reference/reference.fa',
            ".nhr",
            ".nin",
            ".nog",
            ".nsd",
            ".nsi",
            ".nsq"
        )
    conda:
        "../envs/blast.yaml" #because we should not use variables for environments
    log:
        "logs/blastdb.log"
    shell:
        "makeblastdb -in {input} -dbtype nucl -parse_seqids -logfile {log}"


rule mapcontigs:
    input:
        contigs = "denovo_assembly/{sample}/contigs.fasta",
        indexed_ref = multiext(
            'reference/reference.fa',
            ".nhr",
            ".nin",
            ".nog",
            ".nsd",
            ".nsi",
            ".nsq"
        )
    output:
        "blast/contigs/{sample}.tsv"
    log:
        "logs/blast/contigs/{sample}.log"
    threads: 4
    conda:
        "../envs/blast.yaml"
    shell:
        "blastn -query {input.contigs} -db reference/reference.fa -outfmt 6 -out {output} -num_threads {threads} 2>{log}"


rule best_reference:
    input:
        table="blast/contigs/{sample}.tsv",
        reference="reference/reference.fa"
    output:
        "best_references/{sample}.fasta"
    conda:
        "../envs/artificialref.yaml"
    log:
        "logs/best_reference/{sample}.log"
    shell:
        "python scripts/best_reference.py {input.table} {input.reference} {output} 2> {log}"


rule artificialreference:
    input: # contigs for the sample, and best reference for the sample
        best_reference = "best_references/{sample}.fasta", # reference that was most similar to our sample.
        contigs = "denovo_assembly/{sample}/contigs.fasta"
    output:
        artificial_reference = "reference/artificial_reference_{sample}.bam"
    log:
        "logs/artificialreference/{sample}.log"
    threads: 1 #Not possible to assign threads to minimap2
    conda:
        "../envs/artificialref.yaml"
    shell:
        "minimap2 -a {input.best_reference} {input.contigs} | samtools view -b - > {output.artificial_reference} 2> {log}"

rule sort_artificial_reference:
    input:
        "reference/artificial_reference_{sample}.bam"
    output:
        "reference/artificial_reference_{sample}_sorted.bam"
    log:
        "logs/samtools/{sample}_sort.log"
    threads: 4
    conda:
        "../envs/env.yaml"
    shell:
        "samtools sort {input} -o {output} --threads {threads} &> {log}"


rule artificialrefconcensus:
    input:
        best_reference = "best_references/{sample}.fasta",
        bam_sorted = "reference/artificial_reference_{sample}_sorted.bam",
    output:
        consensus = "reference/artificial_reference_{sample}.fa"
    log:
        "logs/bcsf/{sample}.log"
    threads: 1
    conda:
        "../envs/artificialref.yaml"
    shell:
        """
        bcftools mpileup -B -Ou -f {input.best_reference} {input.bam_sorted} | bcftools call -mv -M -Oz -o calls.vcf.gz 2> {log}
        bcftools index calls.vcf.gz -f 2> {log}
        cat {input.best_reference} | bcftools consensus calls.vcf.gz > {output.consensus} 2> {log}
        """
        # TODO; check if it works for you Dana;

        #  [mpileup] 1 samples in 1 input files
        #  [E::bam_plp_push] The input is not sorted (reads out of order)


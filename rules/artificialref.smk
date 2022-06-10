import pandas as pd
configfile : "config/config.yaml"

samples = pd.read_csv(config["samples"],index_col="sample", sep ='\t')

rule makeblastdb:
    input:
        "reference/reference.fa"
    output:
        multiext(
            'reference/references.fa',
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
            'reference/references.fa',
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
        "blastn -query {input.contigs} -db reference/references.fa -outfmt 6 -out {output} -num_threads {threads} 2>{log}"

rule best_reference:
    input:
        table="blast/contigs/{sample}.tsv", reference="reference/references.fa"
    output:
        "best_references/{sample}.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO

        table = pd.read_csv(input[0], sep='\t', header=None)
        table = table[table[0].str.contains('NODE_1')] # get first contig only
        table = table.sort_values([10,2], ascending=[True,False])
        best_ref = table.iloc[0,1]
        with open(input[1]) as ref, open(output,'w') as out:
            record_dict = SeqIO.to_dict(SeqIO.parse(input[1], "fasta"))
            for key in record_dict:
                if best_ref in key:
                    SeqIO.write(record_dict[key], out, "fasta")
                    break







rule readblast:
    input:
        references = config["references"], #fasta file with multiple user provided references, not only prefix.
        contigs = "denovo_assembly/{sample}/contigs.fasta"
    output:
        best_reference = "best_reference_{sample}.fa"
    log:
        "logs/blastreads/{sample}.log"
    threads: 4
    conda:
        "../envs/artificialref.yaml"
    shell:
        "echo hello 2> {log}"

rule artificialreference:
    input: # contigs for the sample, and best reference for the sample
        best_reference = "best_reference_{sample}.fa", # reference that was most similar to our sample.
        contigs = "denovo_assembly/{sample}/contigs.fasta"
    output:
        artificial_reference = "artificial_reference_{sample}.sam"
    log:
        "logs/artificialreference/{sample}.log"
    threads: 4
    conda:
        "../envs/artificialref.yaml"
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
        "logs/bcsf/{sample}.log"
    threads: 4
    conda:
        "../envs/artificialref.yaml"
    shell:
        "bcftools mpileup -B -Ou -f {input.best_reference} {input.sam} | bcftools call -mv -M -Oz -o calls.vcf.gz 2> {log}"
        "bcftools index calls.vcf.gz 2> {log}"
        "cat {input.best_reference} | bcftools consensus calls.vcf.gz > {output.consensus} 2> {log}"
        # TODO; check if it works with 2logs and multiple commands.

rule spades:  # denovo assembly
    input:
        r1 = lambda wildcards: samples.at[wildcards.sample, 'fq1'],
        r2 = lambda wildcards: samples.at[wildcards.sample, 'fq2']
    output:
        file="denovo_assembly/{sample}/contigs.fasta"
    conda:
        "../envs/spades.yaml"
    log:
        "logs/spades/{sample}.log"
    threads: 4
    shell:
        "spades.py -1 {input.r1} -2 {input.r2} -o denovo_assembly/{wildcards.sample} --threads {threads} &>{log}"



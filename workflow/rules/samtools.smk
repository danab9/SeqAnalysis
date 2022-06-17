rule samtobam:
    input:
        "../results/sam/{sample}.sam"
    output:
        "../results/bam/{sample}.bam"
    conda:
       "../envs/env.yaml"
    threads: 4
    log:
        "../results/logs/samtools/{sample}_view.log"
    shell:
      "samtools view -b {input} --threads {threads} > {output} 2> {log}"



rule sort:
    input:
        "../results/bam/{sample}.bam"
    output:
        "../results/bam/sorted/{sample}.bam"
    log:
        "../results/logs/samtools/{sample}_sort.log"
    threads: 4
    conda:
        "../envs/env.yaml"
    shell:
        "samtools sort {input} -o {output} --threads {threads} &> {log}"

rule index:
    input:
        "../results/bam/sorted/{sample}.bam"
    output:
        "../results/bam/sorted/{sample}.bam.bai"
    conda:
        "../envs/env.yaml"
    shell:
        "samtools index -b {input}"

rule mapstats:
    input:
        sorted="../results/bam/sorted/{sample}.bam",
        index="../results/bam/sorted/{sample}.bam.bai"
    output:
        "../results/stats/{sample}.stats"
    threads: 4
    log:
        "../results/log/samtools/{sample}_stats.log"  # added log file
    conda:
        "../envs/env.yaml"
    shell:
        "samtools idxstats {input.sorted} --threads {threads} > {output} 2> {log}"


rule consenus:
    input:
        bam="../results/bam/sorted/{sample}.bam",
        index="../results/bam/sorted/{sample}.bam.bai"
    output:
        "../results/fasta/{sample}.fa"   # "consensus" instead of "fasta/"?
    threads: 1
    log:
        "../results/logs/consenus/{sample}.log"
    conda:
        "../envs/samvar.yaml"
    params:
        q = config["consenus"]["q"]
    shell:
        "samtools mpileup -aa -A -d 0 -Q 0 {input.bam} | ivar consensus -p ../results/fasta/{wildcards.sample} &> {log}"

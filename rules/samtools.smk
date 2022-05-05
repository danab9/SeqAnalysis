configfile : "config/config.yaml"

rule samtobam:
    input:
        "sam/{sample}.sam"
    output:
        "bam/{sample}.bam"
    conda:
       config["env"]
    threads: 4
    log:
        "logs/samtools/{sample}_view.log"
    shell:
      "samtools view -b {input} --threads {threads} > {output} 2> {log}"


rule sort:
    input:
        "bam/{sample}.bam"
    output:
        "bam_sorted/{sample}_sorted.bam"
    log:
        "logs/samtools/{sample}_sort.log"
    threads: 4
    conda:
        config["env"]
    shell:
        "samtools sort {input} -o {output} --threads {threads} &> {log}"

rule index:
    input:
        "bam_sorted/{sample}_sorted.bam"
    output:
        "bam_sorted/{sample}_sorted.bam.bai"
    conda:
        config["env"]
    shell:
        "samtools index -b {input}"

rule mapstats:
    input:
        sorted="bam_sorted/{sample}_sorted.bam",
        index="bam_sorted/{sample}_sorted.bam.bai"
    output:
        "stats/{sample}.stats"
    threads: 4
    conda:
        config["env"]
    shell:
        "samtools idxstats {input.sorted} --threads {threads} > {output}"

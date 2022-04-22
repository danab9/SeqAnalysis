from os import listdir

IDS = [x.split(".")[0] for x in listdir("./sam_tiny")]
sam_dir = "sam_tiny"  # can also get from argv

rule all:
    input: expand("results/{id}.txt", id=IDS)

rule samtobam:
    input:
        sam_dir+"/{dana}.sam"
    output:
        "bam/{dana}.bam"
    shell:
        "samtools view -b {input} > {output}"

rule sort:
    input:
        "bam/{dana}.bam"
    output:
        "bam_sorted/{dana}.bam.sorted"
    shell:
        "samtools sort {input} -o {output}"

rule index:
    input:
        "bam_sorted/{dana}.bam.sorted"
    output:
        "bam_sorted/{dana}.bam.sorted.bai"
    shell:
        "samtools index -b {input}"

rule mapstats:
    input:
        sorted="bam_sorted/{dana}.bam.sorted",
        index="bam_sorted/{dana}.bam.sorted.bai"
    output:
        "results/{dana}.txt"
    shell:
        "samtools idxstats {input.sorted} > {output}"



from os import listdir

sam_dir = "sam_tiny"  # can also get from argv
fastq_dir = "fastq"
refseq = "refseq/reference.fa"
ref_prefix = "refseq/reference"

IDS = [x.split("_")[0]+"_" + x.split("_")[1] for x in listdir(fastq_dir)]

rule all:
    input: expand("results/{id}.rpk", id=IDS)
#     input: expand("sam_tiny/{id}.sam", id=IDS)

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
    threads: 2
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
    shell:
        "bowtie2 -x {ref_prefix} -1 {input.read1} -2 {input.read2} -S {output}"


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

rule rpk:
    input:
         "results/{dana}.txt"
    output:
        "results/{dana}.rpk"
    run:
        with open(input[0]) as file:
            with open(output[0], 'w') as out:
                for line in file:
                    split_line = line.split("\t")
                    length = int(split_line[1])
                    mapped = int(split_line[2])
                    if length != 0:
                        rpk = (mapped/length)/1000
                    else:
                        rpk = 0
                    out.write(line + '\t' + str(rpk))




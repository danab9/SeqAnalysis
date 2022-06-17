def which_fq_1():
    r1=''
    if config['skip_trimming']=="True" and config["decontamination"]!="True":
        r1 = expand(samples.at[sample, 'fq1'], sample=IDS) # (lambda wildcards: samples.at[wildcards.sample, 'fq1'])
    elif config["decontamination"]=='True':
        r1 = "../results/fastq/decontaminated/{sample}_1.fq"
    else:  # config['skip_trimming'] is false
        r1 = "../results/fastq/trimmed/{sample}_1_P.fastq.gz"

    return r1

def which_fq_2():
    r2 = ''
    if config['skip_trimming'] == "True" and config["decontamination"] != "True":
        r2 = expand(samples.at[sample, 'fq2'], sample=IDS) # (lambda wildcards: samples.at[wildcards.sample, 'fq2'])
    elif config["decontamination"] == 'True':
        r2 = "../results/fastq/decontaminated/{sample}_2.fq"
    else:  # config['skip_trimming'] is false
        r2 = "../results/fastq/trimmed/{sample}_2_P.fastq.gz"

    return r2

rule spades:  # denovo assembly
    input:
        # r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1']) if config[
        # "skip_trimming"]=='True' and config["decontamination"] != "True" else "../results/fastq/decontaminated/{sample}_1.fq"
        # if config["decontamination"] == "True" else "../results/fastq/trimmed/{sample}_1_P.fastq.gz",
        # r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2']) if config[
        # "skip_trimming"]=='True' and config["decontamination"] != "True" else "../results/fastq/decontaminated/{sample}_2.fq"
        # if config["decontamination"] == "True" else "../results/fastq/trimmed/{sample}_2_P.fastq.gz"
        r1=which_fq_1(), r2=which_fq_2()
        #r1="../resources/fastq/ERR4082859_1.fastq.gz", r2="../resources/fastq/ERR4082859_2.fastq.gz"

    output:
        "../results/denovo_assembly/{sample}/contigs.fasta"
    conda:
        "../envs/spades.yaml"
    log:
        "../results/logs/spades/{sample}.log"
    threads: 4
    shell:
        "spades.py -1 {input.r1} -2 {input.r2} -o ../results/denovo_assembly/{wildcards.sample} --threads {threads} &>{log}"



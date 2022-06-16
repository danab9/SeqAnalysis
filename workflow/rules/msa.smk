rule msa:
    input:
        sequences = expand("../results/fasta/{sample}.fa", sample=IDS) # added "../results/" before each path
    output:
        alignment = "../results/msa/alignment.fasta"
    log:
        "../results/logs/msa/msa.log"
    threads: 6
    conda:
        "../envs/env.yaml"
    shell:
        "python3 -m augur align --sequences {input.sequences} -o {output.alignment} --threads {threads} &> {log}"

rule tree:
    input:
        alignment = "../results/msa/alignment.fasta"
    output:
        tree = "../results/tree/tree.nwk"
    log:
        "../results/logs/tree/tree.log"
    threads: 6
    params:
        method = config["tree"]["method"]
    conda:
        "../envs/env.yaml"
    shell:
        "augur tree --method {params.method} --alignment {input.alignment} --output {output.tree} --threads {threads} &> {log}"

rule treevisual:
    input:
        tree = "../results/tree/tree.nwk"
    output:
        png = '../results/tree/tree.png'
    conda:
        "../envs/env.yaml"
    threads: 1
    script:
        "../scripts/treevisual.py"

rule variability:
    #  Source: F. Francis, 2015, https://github.com/ffrancis/Multiple-sequence-alignment-Shannon-s-entropy/blob/master/msa_shannon_entropy012915.py !
    input:
        alignment="../results/msa/alignment.fasta"
    params:
        l = config['variability']['l']
    output:
        variability = "../results/variability/variability.txt",
        image = "../results/variability/variability.png" #if image
    log:
        "../results/logs/script/variability.log"
    threads: 1
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/variability.py"

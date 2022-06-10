rule msa:
    input:
        sequences = expand("fasta/{sample}.fa", sample=IDS),
    output:
        alignment = "msa/alignment.fasta"
    log:
        "logs/msa/msa.log"
    threads: 6
    conda:
        "../envs/env.yaml"
    shell:
        "python3 -m augur align --sequences {input.sequences} -o {output.alignment} --threads {threads} &> {log}"

rule tree:
    input:
        alignment = "msa/alignment.fasta"
    output:
        tree = "tree/tree.nwk"
    log:
        "logs/tree/tree.log"
    threads: 6
    params:
        method = config["tree"]["method"]
    conda:
        "../envs/env.yaml"
    shell:
        "augur tree --method {params.method} --alignment {input.alignment} --output {output.tree} --threads {threads} &> {log}"

rule treevisual:
    input:
        tree = "tree/tree.nwk"
    output:
        png = 'tree/tree.png'
    conda:
        "../envs/env.yaml"
    threads: 1
    script:
        "../scripts/treevisual.py"

rule variability:
    #  Source: F. Francis, 2015, https://github.com/ffrancis/Multiple-sequence-alignment-Shannon-s-entropy/blob/master/msa_shannon_entropy012915.py !
    input:
        alignment="msa/alignment.fasta"
    params:
        l = config['variability']['l']
    output:
        variability = "variability/variability.txt",
        image = "variability/variability.png" #if image
    log:
        "logs/script/script.log"
    threads: 1
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/variability.py"

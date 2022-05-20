import pandas as pd
#from ete3 import Tree
configfile : "config/config.yaml"
samples = pd.read_csv(config["samples"],index_col="sample", sep ='\t')
ref_prefix = config['ref']
env = "../envs/env.yaml"

rule msa:
    input:
        sequences = expand("fasta/{sample}.fa", sample=IDS),
        #directory = directory("fasta")
    output:
        alignment = "msa/alignment.fasta"
    log:
        "logs/msa/msa.log"
    threads: 4
    conda:
        env
    shell:
        "python3 -m augur align --sequences {input.sequences} -o {output.alignment} &> {log}"  #Question: is python3 ok?
        #"augur align --sequences {input.directory} -o {output.alignment} &> {log}"



rule tree:
    input:
        alignment = "msa/alignment.fasta"
    output:
        tree = "tree/tree.nwk" #the output does not work yet, but without defining an ouptut, the output is stored as msa/alignment.tree
    log:
        "logs/tree/tree.log"
    threads: 4
    params:
        method = config["tree"]["method"]
    conda:
        env
    shell:
        "augur tree --method {params.method} --alignment {input.alignment} --output {output.tree} &> {log}"     #augur tree --method iqtree --alignment msa/alignment.fasta --output tree/tree.nwk
#
# rule treevisual:
#     input:
#         tree = "tree/tree_raw.nwk"
#     output:
#         png = 'tree.png'
#     run:
#
#     with open('cell_division_newick_cells_0.txt','r') as file:
#         s = file.read()  # read from file
#         t = Tree(s,format=8)  # read Newick format 8
#         print(t)  # print tree as txt
#         t.render('tree.png')
#
# rule treevisual2:
#     # ETE3 http://etetoolkit.org/documentation/ete-view/
#     # https://morpheus.gitlab.io/courses/drawing-cell-genealogies/visualize-trees-using-ete3/
#
#     with open('cell_division_newick_cells_0.txt','r') as file:
#         s = file.read()  # read from file
#         t = Tree(s,format=8)  # read Newick format 8
#         print(t)  # print tree as txt
#         t.render('cell_genealogy.png')


rule variability:
    input:
        alignment="msa/alignment.fasta"
    params:
        l = config['variability']['l']
    output:
        variability = ,
        image = #if image.
    log:
        "logs/script/script.log"
    threads: 4
    conda:
        env
    script:
        "scripts/variability.py"
    #{input.alignment} {output.variability} {params.l} &> {log}"  #Question: is python3 ok?

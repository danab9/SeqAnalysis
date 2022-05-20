from Bio import Phylo
import matplotlib.pyplot as plt

tree_dir = snakemake.input["tree"]
image_dir = snakemake.output["png"]

tree = Phylo.read(tree_dir,"newick")
Phylo.draw(tree)

plt.savefig(image_dir)


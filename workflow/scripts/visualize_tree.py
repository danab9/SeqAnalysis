from biopython import Phylo
import matplotlib.pyplot as plt
#tree_dir = snakemake.input["tree"]
#image_dir = snakemake.output["png"]

tree = Phylo.read("tree/tree.nwk", "newick")
Phylo.draw_ascii(tree)
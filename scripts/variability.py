import pandas as pd

msa = snakemake.input["alignment"]
l = snakemake.params["l"]
image_dir = snakemake.output["image"]
variability_dir = snakemake.output["variability"]

# data = data.sort_values("id")
# data.to_csv(snakemake.output[0], sep="\t")


def entropy(i):
# input: msa column of a poistion
# output: entropy

def entropy_per_i( ):
# input: MSA
# output: entropy per position

function 2:
# input: entropy per position, window width l
# output: averaged entropy over window length (e.g. kernalized)

function 3:
# input: averaged entropy over window length (e.g. kernalized)
# output: visualization.

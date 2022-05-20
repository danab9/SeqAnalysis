from Bio import SeqIO
import numpy as np
import math
import matplotlib.pyplot as plt
msa = snakemake.input["alignment"]
l = snakemake.params["l"]
image_dir = snakemake.output["image"]
variability_dir = snakemake.output["variability"]
# variability_dir = "sldkfja.txt"
# l = 3
# msa = "example.fasta"
# image_dir = "test.png"

alignment_mat = []
with open(msa) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        print(record.id)
        alignment_mat.append([c for c in str(record.seq)]) #vstack

alignment_mat = np.array(alignment_mat)

# data = data.sort_values("id")
# data.to_csv(snakemake.output[0], sep="\t")
#
#


def shannon_entropy(list_input):
    # input: msa column of a poistion
    # output: entropy
    unique_base = set(list_input)  # Get only the unique bases in a column
    # unique_base = unique_base.discard("-")
    M = len(list_input)
    entropy_list = []
    # Number of residues in column
    for base in unique_base:
        n_i = list_input.count(base)  # Number of residues of type i
        P_i = n_i / float(M)  # n_i(Number of residues of type i) / M(Number of residues in column)
        entropy_i = P_i * (math.log(P_i, 2))
        entropy_list.append(entropy_i)
    sh_entropy = -(sum(entropy_list))
    return sh_entropy

def entropy_per_i(alignment_mat):
    # input: MSA
    # output: entropy per position
    return [shannon_entropy(list(alignment_mat[:,i]) ) for i in range(1,alignment_mat.shape[1])]


def window(alignment_mat, l):
    # # input: entropy per position, window width l
    # # output: averaged entropy over window length (e.g. kernalized)
    entropies = entropy_per_i(alignment_mat)
    i =0
    averaged_entopies= []
    while i + l < len(entropies):
        averaged_entopies.append(np.mean(np.array(entropies[i:i+l])))
        i+=1
    return averaged_entopies

window_entropies = window(alignment_mat, l)
with open(variability_dir, "w") as file:
    file.write(str(window_entropies))

plt.plot(window_entropies)
plt.ylabel("entropy")
plt.xlabel("sequence position")
plt.title("Entropy averaged over "+ str(l) + " window length")
plt.savefig(image_dir)
# # input: averaged entropy over window length (e.g. kernalized)
# # output: visualization.

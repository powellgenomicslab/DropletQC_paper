import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import os
import scipy as scipy
import sys as sys

# Define data directory containing loom file and where we will save out other
# objects
data_dir=sys.argv[1]
sample_name=sys.argv[2]
print("Working in directory:")
print(data_dir)
os.chdir(data_dir)

# Import loom file
adata = sc.read_loom(sample_name + ".loom")
# show proportions of spliced/unspliced abundances
scv.utils.show_proportions(adata)

# Get the fraction of UMIs that map to unspliced RNA for each cell and write out
exon_sum = adata.layers['spliced'].sum(axis=1)
intron_sum = adata.layers['unspliced'].sum(axis=1)
nuclear_fraction = intron_sum/(exon_sum + intron_sum)
pd.DataFrame(data=nuclear_fraction,index=adata.obs_names).to_csv("velocyto_nuclear_fraction.csv")

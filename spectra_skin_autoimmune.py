# %%
#import packages
import numpy as np
import json
import scanpy as sc
from collections import OrderedDict
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import pickle

#To track memory
import psutil
import os

#Run command line args
import argparse


#spectra imports
import Spectra as spc
from Spectra import Spectra_util as spc_tl
from Spectra import K_est as kst
from Spectra import default_gene_sets


# %%
#Pull the args
parser = argparse.ArgumentParser(description='Run SPECTRA analysis with subsampled cells')
parser.add_argument('--n_sample_cells', type=int, default=-1, 
                    help='Number of cells to sample for SPECTRA analysis')
parser.add_argument('--n_epochs', type=int, default=2, 
                    help='Number of epochs for SPECTRA training')
parser.add_argument('--bulk_factors', type=bool, default=False, 
                    help='Number of cell factors for SPECTRA training')
parser.add_argument('--n_genes', type=int, default=3000, 
                    help='Number of genes for SPECTRA training')


# Parse the arguments
args = parser.parse_args()

# %% [markdown]
# # Formatting Gene Set Dictionary
# 
# ### Frankensteining SPECTRA's provided dictionary: We used this dictionary to generate the results in the paper: https://doi.org/10.1101/2022.12.20.521311
# 
# ### Options for custom gene_set: We supply the Cytopus knowledge base to construct custom input gene set dictionaries for Spectra. For a tutorial visit the github repository: https://github.com/wallet-maker/cytopus
# 

# %%
#Check memory usage
process = psutil.Process(os.getpid())
print(f"Memory usage: {process.memory_info().rss / 1024**3:.2f} GB")

# %%
#Load skin dataset

skin_data = sc.read_h5ad('/krummellab/data1/DSCoLab/AUTOIPI/papers/skin/data/DB3_54/harmony/sobj_01-13-25.h5ad')
skin_data

# %%
#Add the SPECTRA mapping Emily created
mapping_df = pd.read_csv('~/gabi_rotation/skin_immune_mapping.csv')
spectra_mapping = dict(zip(mapping_df['cell_id'], mapping_df['spectra_annot']))

skin_data.obs['spectra_annotation'] = skin_data.obs.index.map(spectra_mapping)

skin_data.obs['combined_annotation'] = skin_data.obs['spectra_annotation'].fillna(skin_data.obs['cluster_label'])

skin_cats = set(skin_data.obs['combined_annotation'])

# %%
skin_cats

# %%
# load the default gene set dictionary from the Spectra paper:
jpath = '/c4/home/sanjayr/gabi_rotation/Spectra_dict.json'
f = open(jpath, 'r')
annotations = json.loads(f.read())

#Create a new nested dictionary to exactly match our gene set
skin_gene_set = {}

# Add existing valid keys and their values
for key in annotations:
    if key in skin_cats:
        skin_gene_set[key] = annotations[key]

# Add missing keys with empty dictionaries
for key in skin_cats:
    if key not in skin_gene_set:
        skin_gene_set[key] = {}

skin_gene_set['global'] = annotations['global']

# %%
skin_gene_set = spc_tl.check_gene_set_dictionary(
    skin_data,
    skin_gene_set,
    obs_key='combined_annotation',
    global_key='global'
)

# %%
#Check memory usage
process = psutil.Process(os.getpid())
print(f"Memory usage: {process.memory_info().rss / 1024**3:.2f} GB")

# %% [markdown]
# # Create var for highly_variable genes

# %%
#Create boolean mask for top 3000 highly variable genes

n_variable_genes = args.n_genes

most_var_genes = set(skin_data.var['vst.variance'].sort_values(ascending=False)[:n_variable_genes].index)
skin_data.var['highly_variable'] = skin_data.var.index.isin(most_var_genes)

# %% [markdown]
# # Fit SPECTRA model (scaled down)

# %%
#Subsample the skin_data cells

n_sample_cells = args.n_sample_cells
if n_sample_cells == -1 or n_sample_cells > len(skin_data.obs):
    skin_data_sample = skin_data
else:
    skin_data_sample = sc.pp.subsample(skin_data, n_obs=n_sample_cells, copy=True, random_state=42)

skin_gene_set = spc_tl.check_gene_set_dictionary(
    skin_data_sample,
    skin_gene_set,
    obs_key='combined_annotation',
    global_key='global'
)

# %%
# fit the model (We will run this with only 2 epochs to decrease runtime in this tutorial)
eigen_gene_factors = None
if args.bulk_factors:
    with open('factor_estimate.pkl', 'rb') as f:
        eigen_gene_factors = pickle.load(f)
        
n_epochs = args.n_epochs

model = spc.est_spectra(adata=skin_data_sample,
    gene_set_dictionary=skin_gene_set,
    use_highly_variable=True,
    cell_type_key="combined_annotation",
    L = eigen_gene_factors,
    use_weights=True,
    lam=0.1, # varies depending on data and gene sets, try between 0.5 and 0.001
    delta=0.001,
    kappa=None,
    rho=0.001,
    use_cell_types=True,
    n_top_vals=50,
    label_factors=True,
    overlap_threshold=0.2,
    clean_gs = True,
    min_gs_num = 3,
    num_epochs=n_epochs
)


# %%
#Check memory usage
process = psutil.Process(os.getpid())
print(f"Memory usage: {process.memory_info().rss / 1024**3:.2f} GB")

# %%
#save the model uncompressed format
with open(f'models/skin_model_{n_sample_cells}cells_{n_epochs}epochs_bulk{args.bulk_factors}_{n_variable_genes}genes.pickle', 'wb') as f:
    pickle.dump(model, f, pickle.HIGHEST_PROTOCOL)

#save the adata uncompressed format
with open(f'adatas/skin_data_{n_sample_cells}cells_{n_epochs}epochs_bulk{args.bulk_factors}_{n_variable_genes}genes.pickle', 'wb') as f:
    pickle.dump(skin_data_sample, f, pickle.HIGHEST_PROTOCOL)

# At the end of the script, store job metadata
results = {
    "n_sample_cells": n_sample_cells,
    "highly_variable": n_variable_genes,
    "n_epochs": n_epochs,
    "bulk_factors": args.bulk_factors,
}

import json
with open("job_metadata.json", "w") as f:
    json.dump(results, f)

# %% [markdown]
# # LOAD SPECTRA Model

# %%
#and load it like this:
# with open(f'gabi_rotation/models/skin_spectra_model_{n_sample_cells}.pickle', 'rb') as f:
#     model = pickle.load(f)



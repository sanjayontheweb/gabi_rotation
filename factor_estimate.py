# %%
#import packages
import numpy as np
import json
import scanpy as sc
from collections import OrderedDict
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import pickle

#spectra imports
import Spectra as spc
from Spectra import Spectra_util as spc_tl
from Spectra import K_est as kst
from Spectra import default_gene_sets

# %%
skin_data = sc.read_h5ad('/krummellab/data1/DSCoLab/AUTOIPI/papers/skin/data/DB3_54/harmony/sobj_01-13-25.h5ad')

mapping_df = pd.read_csv('~/gabi_rotation/skin_immune_mapping.csv')
spectra_mapping = dict(zip(mapping_df['cell_id'], mapping_df['spectra_annot']))

skin_data.obs['spectra_annotation'] = skin_data.obs.index.map(spectra_mapping)

skin_data.obs['combined_annotation'] = skin_data.obs['spectra_annotation'].fillna(skin_data.obs['cluster_label'])

skin_cats = set(skin_data.obs['combined_annotation'])

# %%
n_variable_genes = 3000

most_var_genes = set(skin_data.var['vst.variance'].sort_values(ascending=False)[:n_variable_genes].index)
skin_data.var['highly_variable'] = skin_data.var.index.isin(most_var_genes)

# %%
L = kst.estimate_L(skin_data, attribute='combined_annotation', highly_variable=True)

# %%
with open('factor_estimate.pkl', 'wb') as f:
    pickle.dump(L, f)

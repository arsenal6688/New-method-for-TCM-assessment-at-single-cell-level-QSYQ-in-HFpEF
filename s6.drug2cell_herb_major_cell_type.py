import scanpy as sc
import drug2cell as d2c
import blitzgsea as blitz
import pandas as pd
import omicverse as ov
import matplotlib.pyplot as plt
import pickle
import blitzgsea as blitz
import os
import numpy as np
import re

output_dir = 'Dir'

adata = ov.read('Dir/celltype_manual.h5ad')
adata

with open("herb-target_dictionary/herb.pkl", "rb") as f:
    chembl_dict_TCM = pickle.load(f)

filtered_dict = chembl_dict_TCM.copy()

sc.pp.normalize_per_cell(adata, counts_per_cell_after = 1e4)
sc.pp.log1p(adata)
adata.raw = adata

hfd_adata = adata[adata.obs['group'] == 'HFD']

d2c.score(hfd_adata, targets=filtered_dict, method="seurat", use_raw=True)

sc.tl.rank_genes_groups(hfd_adata.uns['drug2cell'], method="wilcoxon", groupby="group",groups=['HFD'])

rank_results = sc.get.rank_genes_groups_df(hfd_adata.uns['drug2cell'], group='HFD')

rank_results.to_csv(os.path.join(output_dir, "herb_HFD_rank_results.csv"))

sc.pl.rank_genes_groups_dotplot(hfd_adata.uns['drug2cell'],swap_axes=True, key="rank_genes_groups", dendrogram=False,groupby = 'major_celltype',cmap = 'Blues',categories_order=['Cardiac fibroblast','Endothelial cell','Smooth Muscle Cell','B cell','Neuron','T cell','Pericyte','Macrophage','NKT','Monocyte'],show=False,standard_scale='var')
save_path = os.path.join(output_dir, 'herb_HFD_major_celltype_dotplot.png')
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

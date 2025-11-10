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
from scipy import sparse
import pandas as pd

output_dir = 'Dir'

adata = ov.read('Dir/celltype_manual.h5ad')
adata

with open("ingredient-target_dictionary/HuangQi.pkl", "rb") as f:
    chembl_dict_TCM = pickle.load(f)

num_keys = len(chembl_dict_TCM)
num_keys

target_celltypes = ['Cardiac fibroblast','Endothelial cell',
                    'Macrophage', 'Monocyte']

adata_sub = adata[adata.obs['major_celltype'].isin(target_celltypes)].copy()

d2c.score(adata_sub, targets=chembl_dict_TCM, method="seurat", use_raw=True)

d2c_data = adata_sub.uns['drug2cell']

cell_types = d2c_data.obs['major_celltype'].values
groups = d2c_data.obs['group'].values
gene_names = d2c_data.var_names

hfd_mask = groups == 'HFD'

target_celltypes = ['Cardiac fibroblast', 'Endothelial cell',
                   'Macrophage', 'Monocyte']

results = {}

for cell_type in target_celltypes:
    cell_type_mask = (cell_types == cell_type) & hfd_mask
    if np.sum(cell_type_mask) > 0:  
        cell_type_data = d2c_data.X[cell_type_mask, :]
        if sparse.issparse(cell_type_data):
            fracs = (cell_type_data > 0).mean(axis=0).A1
        else:
            fracs = (cell_type_data > 0).mean(axis=0)
        means = cell_type_data.mean(axis=0)
        if sparse.issparse(means):
            means = means.A1
        results[cell_type] = {
            'fracs': fracs,
            'means': means,
            'cell_count': np.sum(cell_type_mask)
        }

result_dfs = []

for cell_type in target_celltypes:
    if cell_type in results:
        df_temp = pd.DataFrame({
            'gene': gene_names,
            'cell_type': cell_type,
            'frac_hfd': results[cell_type]['fracs'],
            'mean_hfd': results[cell_type]['means'],
            'cell_count': results[cell_type]['cell_count']
        })
        result_dfs.append(df_temp)

final_df = pd.concat(result_dfs, ignore_index=True)
filename = os.path.join(output_dir, f'huangqi_manual_frac_and_mean_count.csv')
final_df.to_csv(filename, index=False)

sorted_results = final_df.sort_values(['cell_type', 'mean_hfd'], ascending=[True, False])

##--
file_path = os.path.join(output_dir, 'huangqi_manual_frac_and_mean_count.csv')
df = pd.read_csv(file_path)

filtered_df = df[
    (df.groupby('gene')['frac_hfd'].transform('max') > 0.25) &
    (df.groupby('gene')['mean_hfd'].transform('max') > 0.5)
]

cell_type_order = ['Cardiac fibroblast', 'Endothelial cell',
                  'Macrophage', 'Monocyte']

max_mean_df = filtered_df.groupby('gene').apply(
    lambda x: x.loc[x['mean_hfd'].idxmax()]
).reset_index(drop=True)

max_mean_df['cell_type'] = pd.Categorical(max_mean_df['cell_type'], 
                                      categories=cell_type_order, 
                                      ordered=True)

max_mean_df_sorted = max_mean_df.sort_values('cell_type')

max_mean_df_sorted.to_csv(os.path.join(output_dir, 'huangqi_manual_frac_and_mean_count_filterd_sorted.csv'), index=False)

top_genes = []
for cell_type in cell_type_order:
    cell_type_mask = max_mean_df_sorted['cell_type'] == cell_type
    top_genes.extend(max_mean_df_sorted[cell_type_mask]['gene'].tolist())

top_genes = list(dict.fromkeys(top_genes))
sc.pl.rank_genes_groups_dotplot(adata_sub.uns['drug2cell'],swap_axes=True, key="rank_genes_groups", dendrogram=False,groupby = 'major_celltype',cmap = 'Blues',standard_scale='var',categories_order=['Cardiac fibroblast', 'Endothelial cell', 'Macrophage','Monocyte'],var_names=top_genes,show=False)
filename = os.path.join(output_dir, f'huangqi_target_celltypes_filtered_components_sorted_dotplot.svg')
plt.savefig(filename, dpi=300, bbox_inches='tight',format='svg')
plt.close()

import omicverse as ov
print(f'omicverse version: {ov.__version__}')
import scanpy as sc
print(f'scanpy version: {sc.__version__}')
ov.ov_plot_set()
import matplotlib.pyplot as plt
import os
import pandas as pd
import pickle

adata = ov.read('Dir/celltype_manual.h5ad')

macrophage_adata = adata[adata.obs['major_celltype'] == 'Macrophage'].copy()
print(f"Number of macrophage: {macrophage_adata.n_obs}")

sc.pp.normalize_total(macrophage_adata)
sc.pp.log1p(macrophage_adata)
sc.pp.highly_variable_genes(macrophage_adata)

sc.pp.pca(macrophage_adata)
sc.pp.neighbors(macrophage_adata, n_neighbors=15, n_pcs=30, use_rep='scaled|original|X_pca')

sc.tl.leiden(macrophage_adata, key_added="macrophage_subtypes", resolution=1)

ov.utils.embedding(macrophage_adata,
                basis='X_mde',
                color=["macrophage_subtypes"],
                title=['Resolution:1'],
                palette=ov.palette()[:],
                frameon='small',show=False
                   )
##--
small_marker_dict={
   'ABCD':['C1QA', 'FCGR1B', 'FCGR1A', 'ADGRE1','TIMD4', 'LYVE1', 'FOLR2', 'IGF1','CCR2', 'ACE', 'PLAC8', 'IL1B','IFIT3', 'IRF7', 'IFIT1B','MMP12']
}

cluster2annotation = {
     '0': 'Timd4+ cells',
     '1': 'Timd4+ cells',
     '2': 'MHC II high cells',
     '3': 'Timd4+ cells',
     '4': 'MHC II high cells',
     '5': 'MHC II high cells',
     '6': 'Timd4+ cells',
     '7': 'MHC II high cells',
     '8': 'CCR2+ cells',
     '9': 'Timd4+ cells',
     '10': 'Timd4+ cells',
     '11': 'interferon-stimulated cells',
     '12': 'MHC II high cells',
     }
macrophage_adata.obs['major_celltype'] = macrophage_adata.obs['macrophage_subtypes'].map(cluster2annotation).astype('category')

high_contrast_palette = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
    '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5'
    ]

macrophage_adata.write("Dir/macrophage_subcelltype_manual.h5ad")

sc.tl.umap(macrophage_adata)

sc.pl.umap(macrophage_adata, color='major_celltype',
          title='UMAP Macrophage',
          palette=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
          show=False, frameon='small')

plt.savefig(
    os.path.join("Dir", "umap_macrophage.svg"),
    dpi=300,
    bbox_inches="tight",
    format='svg'
)
plt.close()

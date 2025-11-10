import omicverse as ov
import scanpy as sc
import os
import matplotlib.pyplot as plt

adata_merged = ov.read('Dir/input_merge_mouse2human.h5ad')

adata_merged=ov.pp.qc(adata_merged,
              tresh={'mito_perc': 0.2, 'nUMIs': 500, 'detected_genes': 200},
              batch_key='batch',
              doublets_method='sccomposite')

sc.pp.calculate_qc_metrics(adata_merged, inplace=True)

min_genes, max_genes = 200, 5000    
adata_merged = adata_merged[
    (adata_merged.obs['n_genes_by_counts'] >= 200) &
    (adata_merged.obs['n_genes_by_counts'] <= 5000)
].copy()

ov.utils.store_layers(adata_merged,layers='counts')

adata_merged=ov.pp.preprocess(adata_merged,mode='shiftlog|pearson',
                       n_HVGs=3000,batch_key='batch')
print(f"Retained cell number: {adata_merged.n_obs}, Gene number: {adata_merged.n_vars}")

adata_merged.raw = adata_merged
adata_merged = adata_merged[:, adata_merged.var.highly_variable_features]

ov.pp.scale(adata_merged)
ov.pp.pca(adata_merged,layer='scaled',n_pcs=50)

adata_merged.obsm["X_mde_pca"] = ov.utils.mde(adata_merged.obsm["scaled|original|X_pca"])

ov.utils.embedding(adata_merged,
                basis='X_mde_pca',frameon='small',
                color=['batch'],show=False)
plt.savefig(
    os.path.join("Dir/", "before_Harmony_for_each_sample.png"),
    bbox_inches="tight",
    dpi=300
)
plt.close()

adata_merged_harmony=ov.single.batch_correction(adata_merged,batch_key='batch',
                                        methods='harmony',n_pcs=50)

adata_merged.obsm["X_mde_harmony"] = ov.utils.mde(adata_merged.obsm["X_harmony"])
ov.utils.embedding(adata_merged,
                basis='X_mde_harmony',frameon='small',
                color=['batch'],show=False)
plt.savefig(
    os.path.join("Dir/", "after_Harmony_for_each_sample.png"),
    bbox_inches="tight",
    dpi=300
)
plt.close()

adata_merged.write("Dir/QC_scale_HVGs_PCA_Harmony.h5ad")

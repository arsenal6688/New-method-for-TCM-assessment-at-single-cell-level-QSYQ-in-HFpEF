import scanpy as sc
import pandas as pd
import numpy as np

adata_merged = sc.AnnData()

adata_GSM7559363 = sc.read_10x_mtx(
    "Dir/GSM7559363",
    var_names="gene_symbols",
    cache=True
)
adata_GSM7559363.obs['batch']='GSM7559363'
adata_GSM7559363.obs['group'] ='HFD'

adata_GSM7559364 = sc.read_10x_mtx(
    "Dir/GSM7559364",
    var_names="gene_symbols",
    cache=True
)
adata_GSM7559364.obs['batch']='GSM7559364'
adata_GSM7559364.obs['group'] ='HFD'

adata_GSM7559365 = sc.read_10x_mtx(
    "Dir/GSM7559365",
    var_names="gene_symbols",
    cache=True
)
adata_GSM7559365.obs['batch']='GSM7559365'
adata_GSM7559365.obs['group'] ='HFD'

adata_GSM7559366 = sc.read_10x_mtx(
    "Dir/GSM7559366",
    var_names="gene_symbols",
    cache=True
)
adata_GSM7559366.obs['batch']='GSM7559366'
adata_GSM7559366.obs['group'] ='WT'

adata_GSM7559367 = sc.read_10x_mtx(
    "Dir/GSM7559367",
    var_names="gene_symbols",
    cache=True
)
adata_GSM7559367.obs['batch']='GSM7559367'
adata_GSM7559367.obs['group'] ='WT'

adata_GSM7559368 = sc.read_10x_mtx(
    "Dir/GSM7559368",
    var_names="gene_symbols",
    cache=True
)
adata_GSM7559368.obs['batch']='GSM7559368'
adata_GSM7559368.obs['group'] ='WT'

adata_merged = sc.concat([adata_GSM7559363, adata_GSM7559364, adata_GSM7559365, adata_GSM7559366, adata_GSM7559367, adata_GSM7559368], axis=0, join='outer', merge='same',
              index_unique="-",keys=['GSM7559363', 'GSM7559364', 'GSM7559365', 'GSM7559366', 'GSM7559367', 'GSM7559368'])
adata_merged.obs['batch'].unique()

total_cells = adata_merged.shape[0]
original_genes = adata_merged.var_names.tolist()
n_original = len(original_genes)
print("Original mouse gene number:", n_original)

homolog_df = pd.read_csv("mouse_human_homologs_dec2021.csv")
mouse_to_human = homolog_df.drop_duplicates("MGI.symbol").set_index("MGI.symbol")["HGNC.symbol"].to_dict()

adata_merged.var["human_gene"] = [mouse_to_human.get(gene, None) for gene in adata_merged.var_names]

valid_genes = ~adata_merged.var["human_gene"].isna()
adata_merged = adata_merged[:, valid_genes].copy()

adata_merged.var_names = adata_merged.var["human_gene"].astype(str) 
del adata_merged.var["human_gene"]
adata_merged.var.index.name = None  

adata_merged = adata_merged[:, ~adata_merged.var_names.duplicated(keep="first")].copy()

sc.pp.filter_genes(adata_merged, min_cells=1)

adata_merged = adata_merged[:, ~adata_merged.var_names.duplicated(keep="first")].copy()

converted_genes = adata_merged.var_names.tolist()
n_converted = len([g for g in converted_genes if g is not None])
print("Converted human gene number:", n_converted)

failed_genes = [g for g,h in zip(original_genes, converted_genes) if h is None]
affected_cells = np.sum(adata_merged[:, failed_genes].X.sum(axis=1) > 0)
print("total raw cell number:", total_cells)
print("affected cell number:", affected_cells)

adata_merged.write("Dir/input_merge_mouse2human.h5ad")

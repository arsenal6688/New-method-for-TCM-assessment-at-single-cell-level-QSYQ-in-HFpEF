import omicverse as ov
print(f'omicverse version: {ov.__version__}')
import scanpy as sc
print(f'scanpy version: {sc.__version__}')
ov.ov_plot_set()
import matplotlib.pyplot as plt
import os
import pandas as pd
import pickle

adata = ov.read('Dir/cluster.h5ad')

adata_raw=adata.raw.to_adata()
adata_raw

##--manual_with_marker_dict
small_marker_dict={
    'Cardiac fibroblast':['PDGFRA','GSN','DCN'],
    'Macrophage':['ADGRE1','C1QA','C1QB'],
    'Smooth Muscle Cell':['TAGLN','ACTA2','MYH11'],
    'Endothelial cell':['CDH5','SOX7','PECAM1'],
    'NKT':['NKG7','GZMA','CCL5'],
    'T cell':['CD3D','CD3E','TRBC1','TRBC2'],
    'B cell':['CD79A','IGHM','IGLC2','IGLL5','IGLL1','IGLC3','IGLC1','IGLC7'],
    'Neuron':['PLP1','KCNA1','KCNA2'],
    'NEUT':['RETNLB','S100A8','S100A9'],
    'Monocyte':['CHIA','PLAC8','CYBB'],
    'Pericyte':['ABCC9','KCNJ8','STEAP4'],
}

smarker_genes_in_data = dict()
for ct, markers in small_marker_dict.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    smarker_genes_in_data[ct] = markers_found

del_markers = list()
for ct, markers in smarker_genes_in_data.items():
    if markers==[]:
        del_markers.append(ct)

for ct in del_markers:
    del smarker_genes_in_data[ct]

smarker_genes_in_data

with open('Dir/smarker_genes_in_data.txt', 'w') as f:
    for cell_type, markers in smarker_genes_in_data.items():
        f.write(f"{cell_type}: {','.join(markers)}\n")

sc.pl.dotplot(
    adata,
    groupby="leiden_res1",
    var_names=smarker_genes_in_data,
    dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    show=False
)
save_path = 'Dir/dotplot_manual_with_marker_dict.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

##--manual_with_marker_diff_genes
adata.uns['log1p']['base']=None
sc.tl.rank_genes_groups(
    adata, groupby="leiden_res1", use_raw=False,
    method="wilcoxon", key_added="dea_leiden_res1"
)

sc.tl.filter_rank_genes_groups(
    adata,
    min_in_group_fraction=0.1,
    max_out_group_fraction=0.2,
    key="dea_leiden_res1",
    key_added="dea_leiden_res1_filtered",
)

sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden_res1", dendrogram=True,
    standard_scale="var", n_genes=3, key="dea_leiden_res1_filtered",
    show=False
)
save_path = 'Dir/dotplot_manual_with_marker_diff_genes.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

clusters = sorted(adata.obs['leiden_res1'].cat.categories)
with open('Dir/all_clusters_top10_markers.txt', 'w') as f:
    for cluster in clusters:
        degs = sc.get.rank_genes_groups_df(
            adata,
            group=str(cluster),
            key='dea_leiden_res1_filtered'
        ).dropna().head(10)
        f.write(f'=== Cluster {cluster} ===\n')
        f.write('\n'.join(degs['names'].astype(str)) + '\n\n')

##--manual_with_marker_database
pathway_dict=ov.utils.geneset_prepare('Dir/CellMarker_Augmented_2021.txt',organism='Human')

adata_aucs=ov.single.pathway_aucell_enrichment(adata_raw,
                                                pathways_dict=pathway_dict,
                                                num_workers=1)
adata_aucs

adata_aucs.obs=adata_raw[adata_aucs.obs.index].obs
adata_aucs.obsm=adata_raw[adata_aucs.obs.index].obsm
adata_aucs.obsp=adata_raw[adata_aucs.obs.index].obsp
adata_aucs

sc.tl.rank_genes_groups(
    adata_aucs, groupby="leiden_res1", use_raw=False,
    method="wilcoxon", key_added="dea_leiden_aucs_res1"
)
sc.pl.rank_genes_groups_dotplot(adata_aucs,groupby='leiden_res1',
                                cmap='RdBu_r',key='dea_leiden_aucs_res1',
                                standard_scale='var',n_genes=3,
                                show=False)
save_path = 'Dir/dotplot_manual_with_marker_database.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

output_path = 'Dir/marker_database_cluster_markers.txt'
with open(output_path, 'w') as f:
    for cluster in adata_aucs.uns['dendrogram_leiden_res1']['categories_ordered']:
        special_cluster = str(cluster)
        degs = sc.get.rank_genes_groups_df(
            adata_aucs,
            group=special_cluster,
            key='dea_leiden_aucs_res1'
        ).dropna()
        f.write(f'{special_cluster}: {"|".join(degs.names[:2].tolist())}\n')

###---Cluster annotation
cluster2annotation = {
     '0': 'Smooth Muscle Cell',
     '1': 'Cardiac fibroblast',
     '2': 'Endothelial cell',
     '3': 'Cardiac fibroblast',
     '4': 'B cell',
     '5': 'Cardiac fibroblast',
     '6': 'Cardiac fibroblast',
     '7': 'Endothelial cell',
     '8': 'Cardiac fibroblast',
     '9': 'Cardiac fibroblast',
     '10': 'Cardiac fibroblast',
     '11': 'Neuron',
     '12': 'Smooth Muscle Cell',
     '13': 'Endothelial cell',
     '14': 'Endothelial cell',
     '15': 'Cardiac fibroblast',
     '16': 'Smooth Muscle Cell',
     '17': 'Endothelial cell',
     '18': 'Cardiac fibroblast',
     '19': 'Smooth Muscle Cell',
     '20': 'Endothelial cell',
     '21': 'Pericyte',
     '22': 'T cell',
     '23': 'Macrophage',
     '24': 'Cardiac fibroblast',
     '25': 'NKT',
     '26': 'Cardiac fibroblast',
     '27': 'T cell',
     '28': 'Monocyte',
     '29': 'Cardiac fibroblast',
     '30': 'Cardiac fibroblast',
     '31': 'Cardiac fibroblast',
}
adata.obs['major_celltype'] = adata.obs['leiden_res1'].map(cluster2annotation).astype('category')

ov.utils.embedding(adata,
                basis='X_mde',
                color=["leiden_res1","major_celltype"],
                title=['Clusters','Major Cell types'],
                palette=ov.palette()[:],wspace=0.55,
                show=False,frameon='small',)

plt.savefig(
    os.path.join("Dir", "manual_major_ct.mde.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

adata.write("Dir/celltype_manual.h5ad")

##--UMAP
sc.tl.umap(adata)

sc.pl.umap(adata, color='major_celltype',
          title='UMAP Major Cell Types',
          palette=ov.palette()[:], 
          show=False, frameon='small')

plt.savefig(
    os.path.join("Dir", "umap_major_celltypes.svg"),
    dpi=300,
    bbox_inches="tight",
    format='svg'
)
plt.close()

import scanpy as sc
import gseapy as gp
from gseapy.plot import barplot, dotplot
import matplotlib.pyplot as plt
import omicverse as ov
import pandas as pd
import numpy as np
import pickle

output_dir = 'Dir'

adata = ov.read('Dir/celltype_manual.h5ad')

##--
groups = ['WT_CD', 'WT_HFD', 'ApoE_KO_CD', 'ApoE_KO_HFD']

sc.tl.rank_genes_groups(
    adata,
    groupby='group',          
    groups=['ApoE_KO_HFD'],           
    reference='WT_CD',           
    method='wilcoxon',        
    corr_method='benjamini-hochberg',  
    tie_correct=True,         
    pts=True                  
)

deg_df = sc.get.rank_genes_groups_df(adata, group='ApoE_KO_HFD')

filter_mask = (
    (deg_df['pvals_adj'] < 0.05) &
    (deg_df['logfoldchanges'].abs() > 1)
)
deg_genes = deg_df[filter_mask]['names'].tolist()

print(f"{len(deg_genes)}")

##--
with open("ingredient-target_dictionary/DanShen.pkl", "rb") as f:
    chembl_dict_TCM = pickle.load(f)

df = pd.read_csv('Dir/danshen_manual_frac_and_mean_count_filterd_sorted.csv')
gene_list = df['gene'].tolist()

all_genes = []

for gene_key in gene_list:
    if gene_key in chembl_dict_TCM:
        all_genes.extend(chembl_dict_TCM[gene_key])

unique_genes = list(set(all_genes))
print(f"{len(unique_genes)}")

overlap_genes = list(set(unique_genes) & set(deg_genes))
print(f"Overlap gene number: {len(overlap_genes)}")

##--
enr = gp.enrichr(
    gene_list=overlap_genes,
    gene_sets=[
        'KEGG_2021_Human',           
        ],
    organism='human'  
)

enr.results.to_csv(f'{output_dir}/KEGG_enrichment.txt', index=False)

sig_pathways = enr.results[
    (enr.results['Adjusted P-value'] < 0.05) &
    (enr.results['Odds Ratio'] > 1) &
    (enr.results['Overlap'].str.split('/').apply(lambda x: int(x[0]) > 5))
].sort_values('Adjusted P-value')

top_20_by_pvalue = sig_pathways.head(20)
sig_pathways = top_20_by_pvalue.sort_values('Odds Ratio', ascending=False)

sig_pathways.to_csv(f'{output_dir}/KEGG_enrichment_sig_top20_pathways.csv', index=False)
print(f"top20: {len(sig_pathways)}")

df = sig_pathways.copy()

df['Gene_ratio'] = df['Overlap'].str.split('/').apply(
    lambda x: int(x[0])/int(x[1])
)
df['-log10(padj)'] = -np.log10(df['Adjusted P-value'])

sig_df = df[
    (df['Adjusted P-value'] < 0.05) &
    (df['Odds Ratio'] > 1)
].sort_values('Odds Ratio').head(20)

sig_df.to_csv(f'{output_dir}/KEGG_enrichment_sig_top20_pathways_filtered.txt', index=False)

plt.figure(figsize=(12,8))
sc = plt.scatter(
    x=sig_df['Odds Ratio'],
    y=sig_df['Term'].str[:50],  
    s=sig_df['Odds Ratio']*5,
    c=sig_df['-log10(padj)'],
    cmap='viridis_r',
    alpha=0.7,
    edgecolor='gray'
)

plt.colorbar(sc, label='-log10(adj.P)')

plt.xlabel('Odds ratio', fontsize=12)

plt.grid(True, linestyle=':', alpha=0.3)
plt.tight_layout()

plt.savefig(f'{output_dir}/KEGG_enrichment.svg', dpi=300, bbox_inches='tight',format='svg')

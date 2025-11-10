import scanpy as sc
import gseapy as gp
from gseapy.plot import barplot, dotplot
import matplotlib.pyplot as plt
import omicverse as ov
import pandas as pd
import numpy as np
from adjustText import adjust_text

output_dir = 'Dir'

adata = ov.read('Dir/celltype_manual.h5ad')

groups = ['WT', 'HFD']

sc.tl.rank_genes_groups(
    adata,
    groupby='group',          
    method='wilcoxon',        
    corr_method='benjamini-hochberg',  
    tie_correct=True,         
    pts=True                  
)

deg_results = sc.get.rank_genes_groups_df(adata, group='HFD')

##--volcano_plot
deg_results['neg_log10_pval'] = -np.log10(deg_results['pvals_adj'])
deg_results['significant_up'] = (deg_results['pvals_adj'] < 0.05) & (deg_results['logfoldchanges'] > 0.25)
deg_results['significant_down'] = (deg_results['pvals_adj'] < 0.05) & (deg_results['logfoldchanges'] < -0.25)

plt.figure(figsize=(12, 8))

non_sig = deg_results[~(deg_results['significant_up'] | deg_results['significant_down'])]
plt.scatter(non_sig['logfoldchanges'], 
           non_sig['neg_log10_pval'],
           c='gray', alpha=0.5, s=8, label='Non-significant')

sig_up = deg_results[deg_results['significant_up']]
plt.scatter(sig_up['logfoldchanges'], 
           sig_up['neg_log10_pval'],
           c='red', alpha=0.7, s=12, label='Up-regulated')

sig_down = deg_results[deg_results['significant_down']]
plt.scatter(sig_down['logfoldchanges'], 
           sig_down['neg_log10_pval'],
           c='blue', alpha=0.7, s=12, label='Down-regulated')

plt.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5, linewidth=1)
plt.axvline(-0.25, color='black', linestyle='--', alpha=0.5, linewidth=1)
plt.axvline(0.25, color='black', linestyle='--', alpha=0.5, linewidth=1)

all_sig_genes = deg_results[(deg_results['significant_up'] | deg_results['significant_down'])]
top_genes = all_sig_genes.nsmallest(10, 'pvals_adj')

plt.xlabel('Log2 Fold Change', fontsize=14)
plt.ylabel('-Log10 Adjusted P-value', fontsize=14)
plt.title('Volcano Plot: DEG', 
          fontsize=16, pad=20)

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

plt.grid(True, alpha=0.3)

n_up = len(sig_up)
n_down = len(sig_down)
n_total = n_up + n_down

textstr = f'Up-regulated: {n_up}\nDown-regulated: {n_down}\nTotal significant: {n_total}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, 
        fontsize=12, verticalalignment='top', bbox=props)

plt.tight_layout()

plt.savefig(f'{output_dir}/volcano_plot_HFD_up_down_colored.png', 
           dpi=300, bbox_inches='tight')

##--Rank plot
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False,show=False)
save_path = 'Dir/sc.tl.rank_genes_groups_top20.svg'
plt.savefig(save_path, dpi=300, bbox_inches='tight',format='svg')
plt.close()

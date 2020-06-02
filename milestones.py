import pandas as pd
import numpy as np

from statsmodels.stats.multitest import multipletests

path_save = '/home/karin/Documents/git/baylor_dicty_paper/try/'
# ************************
# **** Combine all DE (DESeq2) results between neighbouring stages into a single table

# Columns names from DESeq2 result (log2FoldChange, pvalue, and padj) are appended prefix stage1_stage2_
stages = ['no_agg', 'stream', 'lag', 'tag', 'tip', 'slug', 'mhat', 'cul', 'FB']
combined = []
for idx in range(len(stages) - 1):
    stage1 = stages[idx]
    stage2 = stages[idx + 1]
    f = path_save + 'DE_' + stage2 + '_ref_' + stage1 + '.tsv'
    data = pd.read_table(f, index_col=0)[['log2FoldChange', 'pvalue', 'padj']]
    data.columns = [stage1 + '_' + stage2 + '_' + col for col in data.columns]
    combined.append(data)
combined = pd.concat(combined, axis=1, sort=True)

# Adjust FDR across all comparisons. Use only non-nan pvalues.
pvals = []
for col in combined.columns:
    if 'pvalue' in col:
        for row in combined.index:
            pval = combined.at[row, col]
            if not np.isnan(pval):
                pvals.append({'row': row, 'col': col, 'pval': pval})
pvals = pd.DataFrame(pvals)
pvals['FDR_overall'] = multipletests(pvals['pval'], method='fdr_bh')[1]
for row_name, data in pvals.iterrows():
    comparison = data['col'].rstrip('pvalue')
    combined.at[data['row'], comparison + 'FDR_overall'] = data['FDR_overall']

# Save the combined DESeq2 result
combined.to_csv(path_save + 'DE_combined.tsv', sep='\t')

print('Preparing aberrant neighbourhood data.')

import os
import random

from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu, ttest_ind

from helper import *

# *******************
# **** Manage data (load data, specify saving path)
# Path to expression data
path_results = PATH_RESULTS + 'aberrant_neighbourhood/'
if not os.path.exists(path_results):
    os.makedirs(path_results)

# Load expression data
genes = pd.read_csv(PATH_DATA + 'mergedGenes_RPKUM.tsv', sep='\t', index_col=0)
conditions = pd.read_csv(PATH_DATA + 'conditions_mergedGenes.tsv', sep='\t', index_col=None)

# DF with strains and groups
GROUP_DF = []
for strain, group in GROUPS.items():
    GROUP_DF.append({'Strain': strain, 'Group': group})
GROUP_DF = pd.DataFrame(GROUP_DF)

# ************************************************************************
# *** Count in how many replicates per strain the gene is all-zero
genes_rep = merge_genes_conditions(genes=genes, conditions=conditions[['Replicate', 'Measurment']],
                                   matching='Measurment')
genes_rep = split_data(genes_rep, 'Replicate')

genes_zero_count = pd.DataFrame(np.zeros((genes.shape[0], conditions['Strain'].unique().shape[0])),
                                index=genes.index, columns=conditions['Strain'].unique())
for rep, data in genes_rep.items():
    data = data.drop(['Measurment', 'Replicate'], axis=1).T
    strain = conditions[conditions['Replicate'] == rep]['Strain'].values[0]
    data = (data == 0).all(axis=1)
    genes_zero_count[strain] = genes_zero_count[strain] + data
genes_zero_count.to_csv(path_results + 'zero_replicates_count.tsv', sep='\t')

# ***********************************************************
# *** Similarity to closest neighbours of AX4 across strains
# In each strain calculate similarity to the genes that were identified as closest neighbours in AX4
# Neighbours to calculate, this creates 10 neighbours as first neighbours is itself
NEIGHBOURS = 11

# Split data by strain
merged = merge_genes_conditions(genes=genes, conditions=conditions[['Measurment', 'Strain']],
                                matching='Measurment')
splitted = split_data(data=merged, split_by='Strain')
for strain, data in splitted.items():
    splitted[strain] = data.drop(["Strain", 'Measurment'], axis=1).T

# Closest neighbours in AX4
# Exclude genes that are zero in at least one replicate
strain = 'AX4'
print('Calculating neighbours in', strain)
nonzero = genes_zero_count[genes_zero_count[strain] == 0].index
data = splitted[strain]
data = data.loc[nonzero, :]
neighbour_calculator = NeighbourCalculator(genes=data)
neigh_AX4, sims_AX4 = neighbour_calculator.neighbours(n_neighbours=NEIGHBOURS, inverse=False,
                                                      scale=SCALING_SIM, log=LOG_SIM,
                                                      return_neigh_dist=True, remove_self=True)

# Similarity to AX4 neighbours in individual strains on equal-sized vectors.
# Samples are picked at random multiple times, each time calculating similarities on them. Final result is their median.
# Similarities to genes unexpressed in at least one replicate are set to zero. If similarity can not be calculated
# (e.g. whole vector is 0) it is set to 0.
# How many samples to use for each calculation
N_SAMPLES = 10
# How many times recalculate the similarities
N_RESAMPLE = 10
sims_dict = dict()
random.seed(0)
genes_AX4 = set(neigh_AX4.index)
for strain, data in splitted.items():
    print('Calculating neighbours in', strain, 'with resampling')
    nonzero = set(genes_zero_count[genes_zero_count[strain] == 0].index)
    nonzero = list(nonzero & genes_AX4)
    data = data.loc[nonzero, :]
    all_samples = list(data.columns)
    # Repeatedly calculate similarities on sample subset
    similarities_resamples = []
    for i in range(N_RESAMPLE):
        samples = random.sample(all_samples, N_SAMPLES)
        # Keep only data that is not all zero in samples subset
        data_sub = data.loc[:, samples]
        data_sub = data_sub[(data_sub != 0).any(axis=1)]
        # Normalise data for calculation
        data_sub = pd.DataFrame(
            NeighbourCalculator.get_index_query(genes=data_sub, inverse=False, scale='mean0std1', log=True)[0],
            index=data_sub.index, columns=data_sub.columns)
        n_genes = data_sub.shape[0]
        similarities = np.empty((n_genes, NEIGHBOURS - 1))
        similarities[:] = np.nan
        # Calculate similarities for each gene pair, based on the above AX4 results
        for idx_gene in range(n_genes):
            gene = data_sub.index[idx_gene]
            neighbours_WT = neigh_AX4.loc[gene, :].values
            for idx_neigh in range(NEIGHBOURS - 1):
                neigh = neighbours_WT[idx_neigh]
                if neigh in data_sub.index:
                    similarities[idx_gene][idx_neigh] = calc_cosine(data_sub.loc[gene, :], data_sub.loc[neigh, :])
        similarities = pd.DataFrame(similarities, index=data_sub.index, columns=range(NEIGHBOURS - 1))
        # Add genes for which similarities were not calculated, ensures that index is the same for all results
        similarities = similarities.reindex(genes_AX4)
        # Replace non-calculated similarities with 0
        similarities = similarities.replace(np.nan, 0)
        similarities_resamples.append(similarities.values)

    # Calculate median of re-sampled similarities
    similarities_resamples = np.median(np.array(similarities_resamples), axis=0)
    similarities_resamples = pd.DataFrame(similarities_resamples, index=genes_AX4, columns=range(NEIGHBOURS - 1))
    sims_dict[strain] = similarities_resamples

save_pickle(
    path_results + 'sims_AX4basedNeighByStrain_kN' + str(NEIGHBOURS) +
    '_samples' + str(N_SAMPLES) + 'resample' + str(N_RESAMPLE) + '.pkl',
    (neigh_AX4, sims_dict))

# ****************************
# *** Find genes whose similarity to AX4-based neighbourhood is lower in a mutant strain group than in WT

# Make DF with all calculated neighbours of all genes in rows (indexed by gene) and strains in columns.
genes_AX4 = sims_dict['AX4'].index
merged_sims = pd.DataFrame()
for strain, data in sims_dict.items():
    # Make sure that data is in the right format
    if (data.index != genes_AX4).any() or data.isna().any().any():
        raise ValueError('Data index does not match AX4 index or data contains na values')
    merged_sims[strain] = data.values.ravel()
merged_sims.index = np.array([[gene] * sims_dict['AX4'].shape[1] for gene in genes_AX4]).ravel()

# Summarise gene neighbourhood similarity in each strain by median of similarities between gene and neighbours
merged_sims['Gene'] = merged_sims.index
medians = merged_sims.groupby('Gene').median()
merged_sims = merged_sims.drop('Gene', axis=1)
medians.to_csv(path_results + 'simsMedian_AX4basedNeighByStrain_kN' + str(NEIGHBOURS) +
               '_samples' + str(N_SAMPLES) + 'resample' + str(N_RESAMPLE) + '.tsv',sep='\t')

# Pre select genes that are close to their AX4-based neighbourhood in all WT strains
# A threshold for stating that a gene is close to a neighbourhood; based on AX4 data.
threshold = medians['AX4'].quantile(0.2)
print('Threshold for determining if median neighbourhood similarity is high enough in WT:', threshold)
genes_filtered = set(genes_AX4)
for strain in GROUP_DF[GROUP_DF['Group'] == 'WT']['Strain']:
    medians_strain = medians[strain]
    genes_filtered = genes_filtered & set(medians_strain[medians_strain >= threshold].index)
print('N of genes with close AX4-based neighbourhood in WT',len(genes_filtered))

# Compare AX4-based neighbourhood similarity between WT and mutant groups
strains_WT = GROUP_DF[GROUP_DF['Group'] == 'WT']['Strain']
for group in GROUP_DF[GROUP_DF['Group'] != 'WT']['Group'].unique():
    print('Comparing neighbourhood for', group)
    # Strains belonging to a strain group
    strains = GROUP_DF[GROUP_DF['Group'] == group]['Strain']
    group_results = []
    for gene in genes_filtered:
        # Extract similarities to gene neighbours for WT and mutant group
        values_WT = merged_sims.loc[gene, strains_WT].values.ravel()
        values_mutant = merged_sims.loc[gene, strains].values.ravel()

        # Test for mutant group having lower similarities than WT
        result = mannwhitneyu(values_mutant, values_WT, alternative='less')
        p = result[1]
        statistic = result[0]

        m_WT = values_WT.mean()
        m_mutant = values_mutant.mean()
        me_WT = np.median(values_WT)
        me_mutant = np.median(values_mutant)

        group_results.append({'Gene': gene, 'Statistic': statistic, 'p': p,
                              'Mean WT': m_WT, 'Mean mutant': m_mutant, 'Difference mean': m_WT - m_mutant,
                              'Median WT': me_WT, 'Median mutant': me_mutant, 'Difference median': me_WT - me_mutant})
    group_results = pd.DataFrame(group_results)
    # Adjust pvals to FDR
    group_results['FDR'] = multipletests(pvals=group_results['p'], method='fdr_bh', is_sorted=False)[1]
    group_results.to_csv(
        path_results + 'comparisonsSims_' + group + '_AX4basedNeighByStrain_kN' + str(NEIGHBOURS) +
        '_samples' + str(N_SAMPLES) + 'resample' + str(N_RESAMPLE) + '.tsv',
        sep='\t', index=False)

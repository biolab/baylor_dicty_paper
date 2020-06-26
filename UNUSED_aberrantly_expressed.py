import os

from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu, ttest_ind

from helper import *


# *********************************
# *** Helper functions

def quantile_normalise(similarity_means):
    """
    Quantile normalise DF with samples in columns and values in rows, normalise columns to have same distribution.
    Final mapping of sample ranks to rank means is done based on average rank.
    :param similarity_means: DF with samples in columns and values in rows.
    :return: DF of same shape as similarity_means, but with quantile normalised values for each columns.
    """

    # Groupby groups all of the same rank and then averages values for rank
    # Source: https://stackoverflow.com/a/41078786/11521462
    rank_mean = similarity_means.stack().groupby(similarity_means.rank(method='first').stack().astype(int)).mean()
    # Normalise values
    # Find (average) rank and map values to rank-specific values. If between 2 ranks uses their average
    rank_df = similarity_means.rank(method='average')
    quantile_normalised = np.empty(rank_df.shape)
    quantile_normalised[:] = np.nan
    for i in range(rank_df.shape[0]):
        for j in range(rank_df.shape[1]):
            rank = rank_df.iloc[i, j]
            if rank % 1 == 0:
                new_value = rank_mean[rank]
            else:
                rank_low = rank // 1
                rank_high = rank_low + 1
                new_value = (rank_mean[rank_low] + rank_mean[rank_high]) / 2
            quantile_normalised[i, j] = new_value
    normalised = pd.DataFrame(quantile_normalised, index=rank_df.index, columns=rank_df.columns)
    return normalised



# *******************
# **** Manage data (load data, specify saving path)
# Path to expression data
path_results = PATH_RESULTS + 'aberrantly_expressed/'
if not os.path.exists(path_results):
    os.makedirs(path_results)

# Load expression data
genes = pd.read_csv(PATH_DATA + 'mergedGenes_RPKUM.tsv', sep='\t', index_col=0)
conditions = pd.read_csv(PATH_DATA + 'conditions_mergedGenes.tsv', sep='\t', index_col=None)

GROUP_X = {'agg-': 1, 'lag_dis': 2, 'tag_dis': 3, 'tag': 4, 'cud': 5, 'sFB': 6, 'WT': 7, 'prec': 8}

GROUP_DF = []
for strain, group in GROUPS.items():
    GROUP_DF.append({'Strain': strain, 'Group': group, 'X': GROUP_X[group]})
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
sims_dict = dict()
strain = 'AX4'
print('Calculating neighbours in', strain)
nonzero = genes_zero_count[genes_zero_count[strain] == 0].index
data = splitted[strain]
data = data.loc[nonzero, :]
neighbour_calculator = NeighbourCalculator(genes=data)
neigh_AX4, sims_dict[strain] = neighbour_calculator.neighbours(n_neighbours=NEIGHBOURS, inverse=False,
                                                               scale=SCALING_SIM, log=LOG_SIM,
                                                               return_neigh_dist=True, remove_self=True)
# Similarity to AX4 neighbours in other strains
# In individual strains do not calculate similarities for the genes that are all zero in at least one replicate or
# were not used in AX4
for strain, data in splitted.items():
    if strain != 'AX4':
        print('Calculating neighbours in', strain)
        nonzero = set(genes_zero_count[genes_zero_count[strain] == 0].index)
        genes_AX4 = set(neigh_AX4.index)
        nonzero = list(nonzero & genes_AX4)
        data = data.loc[nonzero, :]
        data = pd.DataFrame(
            NeighbourCalculator.get_index_query(genes=data, inverse=False, scale=SCALING_SIM, log=LOG_SIM)[0],
            index=data.index, columns=data.columns)
        n_genes = len(nonzero)
        similarities = np.empty((n_genes, NEIGHBOURS - 1))
        similarities[:] = np.nan
        for idx_gene in range(n_genes):
            gene = nonzero[idx_gene]
            neighbours_AX4 = neigh_AX4.loc[gene, :].values
            for idx_neigh in range(NEIGHBOURS - 1):
                neigh = neighbours_AX4[idx_neigh]
                if neigh in data.index:
                    similarities[idx_gene][idx_neigh] = calc_cosine(data.loc[gene, :], data.loc[neigh, :])
        similarities = pd.DataFrame(similarities, index=nonzero, columns=range(NEIGHBOURS - 1))
        sims_dict[strain] = similarities

save_pickle(
    path_results + 'sims_AX4basedNeighByStrain_kN' + str(NEIGHBOURS) + '.pkl', (neigh_AX4, sims_dict))

# ***********************************************************
# *** Aberrantly expressed genes in strain groups

print('Preparing and normalising neighbours data across strains')
# Make DF with all genes and their neighbours in rows and strains in columns.
# Replace missing similarities (not calculated due to genes having all zero expression) with 0
genes_AX4 = sims_dict['AX4'].index
merged_sims = pd.DataFrame()
for strain, data in sims_dict.items():
    # Reindex individual strain data to have the same gene order across strains.
    data = data.reindex(genes_AX4)
    data = data.replace(np.nan, 0)
    # List all neighbours of all genes in a column, so that neighbours of a gene follow each other
    merged_sims[strain] = data.values.ravel()
merged_sims.index = np.array([[gene] * sims_dict['AX4'].shape[1] for gene in genes_AX4]).ravel()

# Quantile normalise strains data
quantile_normalised = quantile_normalise(similarity_means=merged_sims)
quantile_normalised.to_csv(
    path_results + 'simsQuantileNormalised_AX4basedNeighByStrain_kN' + str(NEIGHBOURS) + '.tsv',
    sep='\t')

# Pre select genes that have high similarity to AX4 based neighbours in WT strains
# Strain specific similarity thresholds obtained for regulons are used
threshold_dict_strain = load_pickle(PATH_RESULTS + 'regulons/' + 'strainThresholds_k2_m0s1log_highest0.7.pkl')
# Retain genes that have median similarity to AX4 neighbours above strain specific threshold in all WT strains
genes_filtered = set(genes_AX4)
print('N genes for which neighbours were calculated in AX4:',len(genes_filtered))
for strain in GROUP_DF[GROUP_DF['Group'] == 'WT']['Strain']:
    threshold = threshold_dict_strain[strain]
    medians = sims_dict[strain].median(axis=1)
    genes_filtered = genes_filtered & set(medians[medians >= threshold].index)
print('N genes with close AX4-based neighbours in all WT strains:',len(genes_filtered))

# For each strain group test if a gene has lower similarity to AX4-based neighbours than in WT strains
strains_WT = GROUP_DF[GROUP_DF['Group'] == 'WT']['Strain']
for group in GROUP_DF[GROUP_DF['Group'] != 'WT']['Group'].unique():
    print('Extracting aberrantly expressed genes for',group)
    strains = GROUP_DF[GROUP_DF['Group'] == group]['Strain']
    group_results = []
    for gene in genes_filtered:
        # Extract similarities to closest neighbours of a gene in WT and mutant group
        values_WT = quantile_normalised.loc[gene, strains_WT].values.ravel()
        values_mutant = quantile_normalised.loc[gene, strains].values.ravel()

        # Test for mutant group similarities being below WT similarities
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
    # Adjust pvals across genes to FDR
    group_results['FDR'] = multipletests(pvals=group_results['p'], method='fdr_bh', is_sorted=False)[1]
    group_results.to_csv(
        path_results + 'comparisonsSims_'+group+'_AX4basedNeigh_u-less_removeZeroRep_scalemean0std1_logTrue_kN11.tsv',
        sep='\t', index=False)

# UNUSED & UNTESTED
# # ***********************************************************
# # *** Aberrantly expressed genes in consecutive strain groups
# # *** Find genes for which neighbours in AX4 do not represent close neighbours in less developed strain groups
#
#
# def compare_gene_scores(quantile_normalised: pd.DataFrame, test: str, alternative: str, group_df: pd.DataFrame,
#                         select_single_comparison: list):
#     """
#     Compare gene scores across strain groups to find genes that have lower/higher score in some set of the groups.
#     It is currently set to work with cosine similarities on the aberrantly expressed genes problem.
#     :param quantile_normalised: DF with strains in columns and genes in rows. For each gene in each strain there is a
#         score.
#     :param select_single_comparison: Order of strain groups used to split groups in two parts.
#         This is a list where strain groups are specified by the value of X column of groups_df.
#         The first element of the list is an ordered list of all strain groups (e.g. ordered by development).
#         The second and third elements are lists of groups that should not be split up into two different groups
#         (can be just the first and the last group from the list of all groups). The splitting is performed between any
#         two neighbouring groups in the specified order so that the reference groups are not split up. Then single
#         possible split is selected based on when the next added strain group would deviate too much from the current
#         high score group.
#     :param test: Test used for testing the difference between groups. 'u' for mannwhitneyu, 't' for t-test.
#     :param alternative: Alternative hypothesis. less (first group has lesser values than the second), greater, two-sided
#     :return: DF with columns: Gene, Comparison (named by the last strain group from the first of the
#         two analysed groups of strain groups), Statistic (test statistic), p (p value), FDR, Mean1, Mean2 (mean
#         of each of two groups from group_splits), Difference (mean2-mean1).
#     """
#     # Split strain groups into two groups. This is done whenever possible, creating multiple possible separations.
#     unsplit1 = select_single_comparison[1]
#     unsplit2 = select_single_comparison[2]
#     select_single_comparison = select_single_comparison[0]
#
#     groups = []
#     g1 = []
#     for group in select_single_comparison:
#         g1.append(group)
#         g2 = [x for x in select_single_comparison if x not in g1]
#         if set(unsplit1).issubset(g1) and set(unsplit2).issubset(g2):
#             groups.append([g1.copy(), g2.copy()])
#
#     # For each gene find best separation into two groups of strain groups and test for difference in gene scores of
#     # strains contained in the two groups
#     results = []
#     for gene in quantile_normalised.index:
#
#         # Find how to split strain groups into two parts
#
#         # Based on the two reference/unsplit groups check which has lower and which higher scores
#         group1 = None
#         group2 = None
#         unsplit_m1 = group_statistic(groups=unsplit1, quantile_normalised=quantile_normalised, gene=gene,
#                                      group_df=group_df, mode='mean')
#         unsplit_m2 = group_statistic(groups=unsplit2, quantile_normalised=quantile_normalised, gene=gene,
#                                      group_df=group_df, mode='mean')
#         if unsplit_m1 > unsplit_m2:
#             order = 1
#         else:
#             order = -1
#         # Select the best group split when next added lower-similarity group starts to deviate too much from the higher
#         # (similarity) set of groups -  too many high's standard deviations away from the high and more than 0.2 away
#         # (expects cosine similarities, this was selected based on aberrantly expressed genes problem).
#         for groups12 in groups[::order]:
#             groups12_ordered = groups12[::order]
#             high_groups = groups12_ordered[0]
#             low_groups = groups12_ordered[1][::order * -1][:-1][::order * -1]
#             group = groups12_ordered[1][::order * -1][-1]
#             if group in unsplit1 or group in unsplit2:
#                 stop = True
#             else:
#                 m = group_statistic(groups=[group], quantile_normalised=quantile_normalised, gene=gene,
#                                     group_df=group_df, mode='mean')
#
#                 m_high = group_statistic(groups=high_groups, quantile_normalised=quantile_normalised, gene=gene,
#                                          group_df=group_df, mode='mean')
#                 m_low = group_statistic(groups=low_groups, quantile_normalised=quantile_normalised,
#                                         gene=gene, group_df=group_df, mode='mean')
#
#                 std_high = group_statistic(groups=high_groups, quantile_normalised=quantile_normalised,
#                                            gene=gene, group_df=group_df, mode='std')
#
#                 group_strains = group_df[group_df['X'] == group]['Strain']
#                 group_values = quantile_normalised.loc[gene, group_strains].values
#
#                 max_diff = 2 * std_high
#                 # Account for very low variability cases in high group - use diff to means instead of 2*std deviations
#                 # from mean if 2*std is smaller than min_diff
#                 min_diff = 0.2
#                 # Change max_diff
#                 if max_diff >= min_diff:
#                     # if any point more than max_diff bellow high mean
#                     stop = (group_values < (m_high - max_diff)).any()
#                 else:
#                     # Diff to means
#                     stop = (m - m_low) <= (m_high - m)
#
#             if stop:
#                 group1 = groups12[0]
#                 group2 = groups12[1]
#                 break
#
#         # Make the two groups based on the above decision of where to split and add name
#         # based on where the split was made
#         compariosn_name = group1[-1]
#         compariosn_name = group_df.query('X == ' + str(compariosn_name)).Group.unique()[0]
#         group_splits = [(group1, group2, compariosn_name)]
#         for comparison in group_splits:
#             # Extract values of both groups
#             strains1 = group_df[group_df['X'].isin(comparison[0])]['Strain']
#             strains2 = group_df[group_df['X'].isin(comparison[1])]['Strain']
#             values1 = quantile_normalised.loc[gene, strains1].values
#             values2 = quantile_normalised.loc[gene, strains2].values
#
#             # Groups summary
#             m1 = values1.mean()
#             m2 = values2.mean()
#             me1 = np.median(values1)
#             me2 = np.median(values2)
#
#             # Test for difference between the two groups
#             if test == 'u':
#                 result = mannwhitneyu(values1, values2, alternative=alternative)
#                 p = result[1]
#                 statistic = result[0]
#             elif test == 't':
#                 result = ttest_ind(values1, values2)
#                 statistic = result[0]
#                 p = result[1]
#                 if alternative == 'less':
#                     if statistic <= 0:
#                         p = p / 2
#                     else:
#                         p = 1 - p / 2
#                 elif alternative == 'greater':
#                     if statistic >= 0:
#                         p = p / 2
#                     else:
#                         p = 1 - p / 2
#             else:
#                 raise ValueError('Test must be t or u.')
#
#             results.append({'Gene': gene, 'Comparison': comparison[2], 'Statistic': statistic, 'p': p,
#                             'Mean1': m1, 'Mean2': m2, 'Difference mean': m2 - m1, 'Difference median': me2 - me1})
#
#     results = pd.DataFrame(results)
#     # Adjust all pvals to FDR
#     results['FDR'] = multipletests(pvals=results['p'], method='fdr_bh', is_sorted=False)[1]
#     return results
#
#
# def group_statistic(groups, quantile_normalised, gene, group_df, mode: str = 'mean'):
#     """
#     Calculates statistic (mean/std) of group of strain groups data (quantile normalised) for a single gene.
#     :param groups: List of strain groups as in X column of GROUP_DF
#     :param quantile_normalised: DF with strains in columns and genes in rows. For each gene in each strain there is a
#         score (e.g. average similarity to neighbours).
#     :param gene: Gene ID for which to calculate the statistic
#     :param mode: mean or std or median
#     """
#     strains = group_df[group_df['X'].isin(groups)]['Strain']
#     values = quantile_normalised.loc[gene, strains].values
#     if mode == 'mean':
#         return values.mean()
#     elif mode == 'std':
#         return values.std()
#     elif mode == 'median':
#         return values.median()
#
# # Calculate average similarity to the AX4 based neighbours
# # For some genes pairs similarities could not be calculated - set these similarities to zero
# # This is also done for genes for which no neighbours were calculated
# similarity_means = pd.DataFrame(index=genes.index, columns=sims_dict.keys())
# for strain, data in sims_dict.items():
#     data = data.replace(np.nan, 0)
#     means = data.reindex(similarity_means.index).mean(axis=1)
#     means = means.replace(np.nan, 0)
#     similarity_means[strain] = means
#
# # Quantile normalise average similarities to AX4 based neighbours in each strain
# quantile_normalised = quantile_normalise(similarity_means=similarity_means)
#
# quantile_normalised.to_csv(
#     path_results + 'quantileNormAvgSims_AX4basedNeighByStrain_kN' + str(NEIGHBOURS) + '.tsv', sep='\t')
#
# # Pre select genes that do not have low similarity to AX4 neighbours in developing strains (WT, Prec)
# genes_filtered = set(similarity_means.index)
# for strain in GROUP_DF[GROUP_DF['Group'].isin(['prec', 'WT'])]['Strain']:
#     threshold = np.quantile(similarity_means[strain], 0.3)
#     genes_filtered = genes_filtered & set(similarity_means[similarity_means[strain] >= threshold].index)
#
# # *** Find genes that have low similarity to AX4 neighbours in less developed and high in more developed strains
#
# test = 'u'
# alternative = 'less'
#
# # Compare less developed strain groups to more developed.
# # Do not use sFB. Treat both disaggregation groups as a single disaggregation group.
# group_df = GROUP_DF.copy()
# group_df.loc[group_df['X'] == 2, 'X'] = 3
# group_df.loc[group_df['X'] == 3, 'Group'] = 'dis'
# select_single_comparison = [[1, 3, 4, 5, 7, 8], [1], [7, 8]]
#
# results_WT = compare_gene_scores(quantile_normalised=quantile_normalised.loc[genes_filtered, :], test=test,
#                                  alternative=alternative, select_single_comparison=select_single_comparison,
#                                  group_df=group_df)
# results_WT.to_csv(
#     path_results + 'comparisonsAvgSims_AX4basedNeigh_' + test + '-' + alternative + '.tsv',
#     sep='\t', index=False)

print('Preparing regulons data.')

import glob
import os

from helper import *


# *******************
# *** Manage data (load data, specify saving path)
# Path to expression data
path_results = PATH_RESULTS + 'regulons/'
if not os.path.exists(path_results):
    os.makedirs(path_results)

# Load expression data
genes = pd.read_csv(PATH_DATA + 'mergedGenes_RPKUM.tsv', sep='\t', index_col=0)
conditions = pd.read_csv(PATH_DATA + 'conditions_mergedGenes.tsv', sep='\t', index_col=None)

# *************************
# *** Similarity thresholds for each replicate
SPLITBY = 'Strain'
merged = merge_genes_conditions(genes=genes, conditions=conditions[['Measurment', SPLITBY]],
                                matching='Measurment')
splitted = split_data(data=merged, split_by=SPLITBY)
for rep, data in splitted.items():
    splitted[rep] = data.drop([SPLITBY, 'Measurment'], axis=1).T

strains = conditions['Strain'].unique()

# Calculate similarities to the closest neighbour
sims_dict = dict()
for strain in strains:
    print('Calculating closest 2 neighbours for', strain)
    data = splitted[strain]
    neighbour_calculator = NeighbourCalculator(genes=data)
    neigh, sims_dict[strain] = neighbour_calculator.neighbours(n_neighbours=2, inverse=False, scale=SCALING_SIM,
                                                               log=LOG_SIM, return_neigh_dist=True, remove_self=True)

# Obtain thresholds from similarities to the closest neighbour across genes (30th percentile) of each strain
THRESHOLD_QUANTILE = 0.3
threshold_data = {}
for strain in strains:
    # Only 2 neighbours were calculated and one of them (self) was removed, so only one neighbour is left
    sims = sims_dict[strain][0]
    # Thresholds were rounded for display (not shown) and subsequent analysis
    threshold = np.round(np.quantile(sims, THRESHOLD_QUANTILE), 3)
    threshold_data[strain] = threshold

save_pickle(path_results + 'strainThresholds_k2_m0s1log_highest' + str(1 - THRESHOLD_QUANTILE) + '.pkl', threshold_data)

# ***********************
# *** Find close gene neighbours in each strain

# Find samples of each strain
batches = list(conditions['Strain'])
batches = np.array(batches)

# Number of neighbours to calculate for each gene.
KNN = 300
# Calculate neighbours
for batch in set(batches):
    print('Calculating neighbours for', batch)
    genes_sub = genes.T[batches == batch].T
    neighbour_calculator = NeighbourCalculator(genes_sub)
    result = neighbour_calculator.neighbours(KNN, inverse=False, scale=SCALING_SIM, log=LOG_SIM)
    save_pickle(path_results + 'neighbours_kN' + str(KNN) + '_mean0std1_log_' + batch + '.pkl', result)

# **************************
# *** Combine close neighbour results across strains

# Extract genes from strains and merge into a single matrix
files = [f for f in glob.glob(path_results + 'neighbours_kN' + str(KNN) + '_mean0std1_log_' + "*.pkl")]
n_genes = genes.shape[0]
genes_dict = dict(zip(genes.index, range(n_genes)))
merged_results = np.zeros((n_genes, n_genes))
for f in files:
    strain = f.split('_')[-1].replace('.pkl', '')
    result = load_pickle(f)
    similarity_threshold = threshold_data[strain]
    print('Adding data of', strain, ' using similarity threshold', similarity_threshold)
    for pair, similarity in result.items():
        # If similarity of a gene pair is at least equal to the similarity threshold then count the two genes
        # as co-expressed in a strain
        if similarity >= similarity_threshold:
            gene1 = genes_dict[pair[0]]
            gene2 = genes_dict[pair[1]]
            merged_results[gene1, gene2] += 1
            merged_results[gene2, gene1] += 1
merged_results = pd.DataFrame(merged_results, index=genes.index, columns=genes.index)

merged_results.to_csv(path_results + 'strainMergedNeighbours_kN' + str(KNN) + '_mean0std1_log_highest' +
                      str(1 - THRESHOLD_QUANTILE) + '.tsv', sep='\t')

# ***************************
# ***** Extract genes with close neighbours in strains in which the gene is expressed

# Expression threshold for assuming that a gene is expressed in a strain:
# 10% of 99th percentile of expression across samples
ratio_max = 0.99
proportion = 0.1
threshold = genes.quantile(q=ratio_max, axis=1) * proportion

# For each gene count how many strains are assumed to have the gene expressed, based on the above threshold
SPLITBY = 'Strain'
merged = merge_genes_conditions(genes=genes, conditions=conditions[['Measurment', SPLITBY]],
                                matching='Measurment')
splitted = split_data(data=merged, split_by=SPLITBY)
for rep, data in splitted.items():
    splitted[rep] = data.drop([SPLITBY, 'Measurment'], axis=1).T

strain_expressed = pd.DataFrame()
for strain, data in splitted.items():
    strain_expressed[strain] = (data.T >= threshold).any()
n_strains = strain_expressed.sum(axis=1)

# Has a close neighbour in at least one strain
min_strains = 1
n_strains[n_strains < min_strains] = min_strains
# Do not require that a gene has a neighbour in all strains where it is expressed (max possible 21), but rather
# set limit at 18 to reduce false negatives. On the other hand, genes termed as expressed in less strains (<21)
# may actually have neighbours in more than just these strains.
max_strains = 18
n_strains[n_strains > max_strains] = max_strains

# For each gene: maximal number of strains in which a close neighbour is present
genemax = merged_results.max()

# Remove genes whose neighbours are not present in at least N strains, using gene specific N
remove_genes = genemax.loc[genemax < n_strains].index
merged_results_filtered = merged_results.drop(index=remove_genes, columns=remove_genes)

merged_results_filtered.to_csv(path_results + 'regulonGenes_kN' + str(KNN) + '_mean0std1_log_highest' +
                               str(1 - THRESHOLD_QUANTILE) +
                               '_minExpressed' + str(ratio_max) + str(proportion) +
                               '_strainFilterMin' + str(min_strains) + 'Max' + str(max_strains) + '.tsv', sep='\t')

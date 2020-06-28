import os

from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu, ttest_ind

from helper import *


# *******************
# **** Manage data (load data, specify saving path)
# Path to expression data
path_results = PATH_RESULTS + 'aberrantly_expressed/'
if not os.path.exists(path_results):
    os.makedirs(path_results)

# Load expression data
genes = pd.read_csv(PATH_DATA + 'mergedGenes_RPKUM.tsv', sep='\t', index_col=0)
conditions = pd.read_csv(PATH_DATA + 'conditions_mergedGenes.tsv', sep='\t', index_col=None)

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


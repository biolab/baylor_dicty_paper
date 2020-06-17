print('Preparing regulons data.')

import warnings
import glob
import os

import multiprocessing
import numpy as np
import pandas as pd

from pynndescent import NNDescent
import sklearn.preprocessing as pp

from helper import merge_genes_conditions, split_data, save_pickle, load_pickle, PATH_RESULTS, PATH_DATA

# ***************
# *** Helper functions
SCALING = 'mean0std1'
LOG = True

THREADS = multiprocessing.cpu_count() - 1
if THREADS > 30:
    THREADS = 30


class NeighbourCalculator:
    """
    Obtains best neighbours of genes based on their expression profile.
    """

    MINMAX = 'minmax'
    MEANSTD = 'mean0std1'
    NONE = 'none'
    SCALES = [MINMAX, MEANSTD, NONE]

    def __init__(self, genes: pd.DataFrame, remove_zero: bool = True):
        """
        :param genes: Data frame of genes in rows and samples in columns. Index is treated as gene names
        :param remove_zero: Remove genes that have all expression values 0.
        """
        NeighbourCalculator.check_numeric(genes)
        if remove_zero:
            genes = genes[(genes != 0).any(axis=1)]
        self._genes = genes

    @staticmethod
    def check_numeric(expression: pd.DataFrame):
        """
        Is content of data frame numeric (pandas float/int).
        If data is not numeric raises exception.
        :param expression: data frame to be analysed
        """
        if not ((expression.dtypes == 'float64').all() or (expression.dtypes == 'int64').all()):
            raise ValueError('Expression data is not numeric.')

    def neighbours(self, n_neighbours: int, inverse: bool, scale: str = SCALING, log: bool = LOG,
                   return_neigh_dist: bool = False, genes_query_names: list = None, remove_self: bool = False):
        """
        Calculates neighbours of genes based on expression profile across samples (columns) based on cosine distances.
        Wrapper for calculate_neighbours.
        :param n_neighbours: Number of neighbours to obtain for each gene. This will include self for non-inverse.
        :param inverse: Calculate most similar neighbours (False) or neighbours with inverse profile (True).
        :param scale: Scale expression by gene with 'minmax' (min=0, max=1) or 'mean0std1' (mean=0, std=1) or 'none'.
        :param log: Should expression data be log2(data+pseudocount) transformed before scaling.
        :param return_neigh_dist: Return tuple with nearest neighbour matrix and similarity matrix data frames,
            as returned by pynndescent, but with similarity matrix converted to similarities and with added gene
            names for the index.
        :param genes_query_names: Use only the specified genes as query (out of all genes).
        :param remove_self: Used only if return_neigh_dist is true. Whether to remove sample from its closest
            neighbours or not. If return_neigh_dist is False this is done automatically. This also removes the last
            column of neighbours if self is not present - thus it should not be used with inverse,
            as self will not be present.
        :return: Dict with keys being gene pair names tuple (smaller name by alphabet is the first tuple value) and
            values representing cosine similarity. Or see return_neigh_dist.
        """
        if scale not in self.SCALES:
            raise ValueError('Scale must be:', self.SCALES)

        genes = self._genes

        if genes_query_names is not None:
            genes_query = genes.loc[genes_query_names, :]
        else:
            genes_query = None
        return NeighbourCalculator.calculate_neighbours(genes=genes, n_neighbours=n_neighbours, inverse=inverse,
                                                        scale=scale,
                                                        log=log, return_neigh_sim=return_neigh_dist,
                                                        genes_query_data=genes_query, remove_self=remove_self)

    @staticmethod
    def calculate_neighbours(genes, n_neighbours: int, inverse: bool, scale: str, log: bool,
                             description: str = '', return_neigh_sim: bool = False,
                             genes_query_data: pd.DataFrame = None, remove_self: bool = False):
        """
        Calculate neighbours of genes based on cosine distance.
        :param genes: Data frame as in class init, gene names (rows) should match the one in init.
        :param n_neighbours: Number of neighbours to obtain for each gene. This will include self for non-inverse.
        :param inverse: Calculate most similar neighbours (False) or neighbours with inverse profile (True).
        :param scale: Scale expression by gene with 'minmax' (min=0, max=1) or 'mean0std1' (mean=0, std=1) or 'none'.
        :param log: Should expression data be log2(data+pseudocount) transformed before scaling.
        :param description: If an error occurs while making KNN index report this description with the error.
        :param return_neigh_sim: Return tuple with nearest neighbour matrix and similarity matrix data frames,
            as returned by pynndescent, but with distance matrix converted to similarities and with added gene
            names for the index.
        :param genes_query_data: Use this as query. If None use genes.
        :param remove_self: Used only if return_neigh_dist is true. Whether to remove sample from its closest
            neighbours or not. If return_neigh_dist is False this is done automatically. This also removes the last
            column of neighbours if self is not present - thus it should not be used with inverse,
            as self will not be present.
        :return: Dict with keys being gene pair names tuple (smaller name by alphabet is the first tuple value) and
            values representing cosine similarity. Or see return_neigh_dist.
        """
        genes_index, genes_query = NeighbourCalculator.get_index_query(genes=genes, inverse=inverse, scale=scale,
                                                                       log=log,
                                                                       genes_query_data=genes_query_data)
        # Random state was not set during the analysis in the paper so the obtained results might differ slightly
        try:
            index = NNDescent(genes_index, n_jobs=THREADS, metric='cosine', random_state=0)
        except ValueError:
            try:
                index = NNDescent(genes_index, tree_init=False, n_jobs=THREADS, random_state=0)
                warnings.warn(
                    'Dataset ' + description + ' index computed without tree initialisation',
                    Warning)
            except ValueError:
                raise ValueError('Dataset ' + description + ' can not be processed by pydescent')
        neighbours, distances = index.query(genes_query.tolist(), k=n_neighbours)

        if genes_query_data is None:
            genes_query_data = genes
        if return_neigh_sim:
            neighbours = NeighbourCalculator.parse_neighbours_matrix(neighbours=neighbours,
                                                                     genes_query=genes_query_data,
                                                                     genes_idx=genes)
            similarities = pd.DataFrame(NeighbourCalculator.parse_distances_matrix(distances),
                                        index=genes_query_data.index)
            if remove_self:
                neighbours, similarities = NeighbourCalculator.remove_self_pynn_matrix(neighbours=neighbours,
                                                                                       similarities=similarities)
            return neighbours, similarities
        else:
            return NeighbourCalculator.parse_neighbours(neighbours=neighbours, distances=distances,
                                                        genes_query=genes_query_data, genes_idx=genes)

    @classmethod
    def get_index_query(cls, genes: pd.DataFrame, inverse: bool, scale: str = SCALING, log: bool = LOG,
                        genes_query_data: pd.DataFrame = None) -> tuple:
        """
        Get genes data scaled to be index or query for neighbour search.
        :param genes: Gene data for index and query.
        :param inverse: If True inverse query to compute neighbours with opposite profile.
        :param scale: Scale expression by gene with 'minmax' (min=0, max=1) or 'mean0std1' (mean=0, std=1) or 'none'.
        :param log: Should expression data be log2(data+pseudocount) transformed.
        :param genes_query_data: Genes data for query, if None uses genes.
        :return: genes for index (1st element) and genes for query (2nd element)
        """
        if scale not in [cls.MINMAX, cls.MEANSTD, cls.NONE]:
            raise ValueError('This scaling parameter is unknown')
        if log:
            genes = np.log2(genes + 1)
            if genes_query_data is not None:
                genes_query_data = np.log2(genes_query_data + 1)
        if inverse:
            if genes_query_data is None:
                genes_query_data = genes
            genes_query = genes_query_data * -1
            genes_index = genes
            if scale == cls.MINMAX:
                genes_query = cls.minmax_scale(genes_query)
                genes_index = cls.minmax_scale(genes_index)
            elif scale == cls.MEANSTD:
                genes_query = cls.meanstd_scale(genes_query)
                genes_index = cls.meanstd_scale(genes_index)
            elif scale == cls.NONE:
                genes_query = genes_query.values
                genes_index = genes_index.values
        else:
            if scale == cls.MINMAX:
                genes = cls.minmax_scale(genes)
                if genes_query_data is not None:
                    genes_query_data = cls.minmax_scale(genes_query_data)
            elif scale == cls.MEANSTD:
                genes = cls.meanstd_scale(genes)
                if genes_query_data is not None:
                    genes_query_data = cls.meanstd_scale(genes_query_data)
            elif scale == cls.NONE:
                genes = genes.values
                if genes_query_data is not None:
                    genes_query_data = genes_query_data.values
            genes_index = genes
            if genes_query_data is None:
                genes_query_data = genes
            genes_query = genes_query_data
        return genes_index, genes_query

    @staticmethod
    def minmax_scale(genes: pd.DataFrame) -> np.ndarray:
        """
        Scale each row from 0 to 1.
        :param genes: data
        :return: scaled data
        """
        return pp.minmax_scale(genes, axis=1)

    @staticmethod
    def meanstd_scale(genes: pd.DataFrame) -> np.ndarray:
        """
        Scale each row to mean 0 and std 1.
        :param genes: data
        :return: scaled data
        """
        return pp.scale(genes, axis=1)

    @staticmethod
    def parse_neighbours(neighbours: np.ndarray, distances: np.ndarray, genes_idx: pd.DataFrame,
                         genes_query: pd.DataFrame) -> dict:
        """
        Transform neighbours and distances data frames into dictionary with neighbours (gene pair names
        sorted alphabetically) as keys and values being similarities.
        If pair of neighbours is given more than once (e.g. twice) the added similarity is an average of
        the already present and the new similarity.
        The neighbours are named based on genes data (index) stored in NeighbourCalculator instance.
        :param neighbours: Array of shape: genes*neighbours, where for each gene there are specified neighbours
            (as indices).
        :param distances: Array of shape: genes*neighbours, where for each gene there are specified distances
            for each neighbour, in the same order as in neighbours DF.
        :param genes_query: DataFrame used to name genes in neighbours & distances rows (the query of PyNNDescent),
            gene names in rows.
        :param genes_idx: DataFrame used to name genes in neighbours & distances table entries (the PyNNDescent index),
            gene names in rows.
        :return: Dict with gene names as tuple keys (smaller by alphabet is first tuple value) and
            values representing similarity: 1-distance
        """
        parsed = dict()
        for gene in range(distances.shape[0]):
            for neighbour in range(distances.shape[1]):
                distance = distances[gene, neighbour]
                similarity = NeighbourCalculator.cosine_dist_to_sim(distance)
                gene2 = neighbours[gene, neighbour]
                gene_name1 = genes_query.index[gene]
                gene_name2 = genes_idx.index[gene2]
                if gene_name1 != gene_name2:
                    if gene_name2 > gene_name1:
                        add_name1 = gene_name1
                        add_name2 = gene_name2
                    else:
                        add_name1 = gene_name2
                        add_name2 = gene_name1
                    if (add_name1, add_name2) in parsed.keys():
                        # Can do averaging directly with existing value as there will not be more than 2 pairs
                        # with the same genes (eg. both possible positions: a,b and b,a for index,query)
                        # Similarities may be different when inverse is used with minmax
                        parsed[(add_name1, add_name2)] = (parsed[(add_name1, add_name2)] + similarity) / 2
                    else:
                        parsed[(add_name1, add_name2)] = similarity
        return parsed

    @staticmethod
    def parse_neighbours_matrix(neighbours: np.ndarray, genes_idx: pd.DataFrame,
                                genes_query: pd.DataFrame) -> pd.DataFrame:
        """
        Names pynndescent neighbours table (values and rows).
        :param neighbours: As returned by pynndescent
        :param genes_idx: Data frame with rownames matching neighbours values (used for pynndescent index)
        :param genes_query: Data frame with rownames matching neighbours rows (used for pynndescent query)
        :return: Named neighbours table.
        """
        parsed = pd.DataFrame(columns=range(neighbours.shape[1]))
        for gene1 in range(neighbours.shape[0]):
            for gene2_col in range(neighbours.shape[1]):
                gene2 = neighbours[gene1, gene2_col]
                gene_name1 = genes_query.index[gene1]
                gene_name2 = genes_idx.index[gene2]
                parsed.loc[gene_name1, gene2_col] = gene_name2
        return parsed

    @staticmethod
    def parse_distances_matrix(distances: np.ndarray):
        """
        Transform cosine distances to similarities
        :param distances: pynndescent cosine distance matrix
        """
        return NeighbourCalculator.cosine_dist_to_sim(distances)

    @staticmethod
    def remove_self_pynn_matrix(neighbours: pd.DataFrame, similarities: pd.DataFrame):
        """
        Parse similarities and neighbours data frames. Remove entries of each gene that represent itself being
        a close neighbour. If self is not in neghbours then it removes the last element from the neighbours and
        similarities, so that all genes have same N of neighbours.
        :param similarities: Similarities matrix from neighbours function (row names with gene names)
        :param neighbours: Neighbours matrix from neighbours function (row names and values named with genes)
        """
        similarities_parsed = pd.DataFrame(index=similarities.index, columns=similarities.columns[:-1])
        neighbours_parsed = pd.DataFrame(index=neighbours.index, columns=neighbours.columns[:-1])
        for gene in neighbours.index:
            neigh_row = neighbours.loc[gene, :]
            sims_row = similarities.loc[gene, :]
            if gene not in neigh_row.values:
                similarities_parsed.loc[gene, :] = sims_row.iloc[:-1]
                neighbours_parsed.loc[gene, :] = neigh_row.iloc[:-1]
            else:
                self_idx = neigh_row[neigh_row == gene].index[0]
                similarities_parsed.loc[gene, :] = sims_row.drop(similarities.columns[self_idx]).values
                neighbours_parsed.loc[gene, :] = neigh_row.drop(neighbours.columns[self_idx]).values
        return neighbours_parsed, similarities_parsed

    @staticmethod
    def cosine_dist_to_sim(dist):
        """
        Works on cosine and Pearson correlation distances.
        An exception is raised if the number is  different from 0 or 2 after rounding to 4 decimal places.
        :param dist: Float or np.ndaray
        """
        # Some distances are slightly below 0 or above 2 due to numerical precision.
        if isinstance(dist, float):
            if round(dist, 4) < 0:
                raise ValueError('Odd cosine distance below 0.')
            if round(dist, 4) > 2:
                raise ValueError('Odd cosine distance above 2.')
        elif isinstance(dist, np.ndarray):
            if np.around(dist, 4).any() < 0:
                raise ValueError('Odd cosine distance below 0.')
            if np.around(dist, 4).any() > 2:
                raise ValueError('Odd cosine distance above 2.')
        else:
            raise ValueError('The parameter dist must be float or np.ndarray.')
        return 1 - dist


# *******************
# **** Manage data (load data, specify saving path)
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
    neigh, sims_dict[strain] = neighbour_calculator.neighbours(n_neighbours=2, inverse=False, scale=SCALING,
                                                               log=LOG, return_neigh_dist=True, remove_self=True)

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
    result = neighbour_calculator.neighbours(KNN, inverse=False, scale=SCALING, log=LOG)
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

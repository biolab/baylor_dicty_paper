import warnings
import glob
import pickle

import numpy as np
import pandas as pd

from pynndescent import NNDescent
import sklearn.preprocessing as pp

# ***************
# *** Helper functions
# TODO check helper functions (arguments, descriptions)
SCALING = 'mean0std1'
LOG = True


def save_pickle(file: str, object):
    """
    Pickle an object into a file.
    :param file: Absolute file path.
    :param object: Object to pickle.
    """
    f = open(file, 'wb')
    pickle.dump(object, f)
    f.close()


def load_pickle(file: str):
    """
    Read a pickled object.
    :param file: Absolute file name.
    :return: Object.
    """
    pkl_file = open(file, 'rb')
    result = pickle.load(pkl_file)
    pkl_file.close()
    return result


class NeighbourCalculator:
    """
    Obtains best neighbours of genes based on their expression profile.
    This can be done for all conditions at once or by condition groups, later on merging the results into
    single neighbour data.
    """

    MINMAX = 'minmax'
    MEANSTD = 'mean0std1'
    NONE = 'none'
    SCALES = [MINMAX, MEANSTD, NONE]

    def __init__(self, genes: pd.DataFrame, remove_zero: bool = True, conditions: pd.DataFrame = None,
                 conditions_names_column=None):
        """
        :param genes: Data frame of genes in rows and conditions in columns. Index is treated as gene names
        :param remove_zero: Remove genes that have all expression values 0.
            If batches are latter specified there may be all 0 rows within individual batches.
        :param conditions: data frame with conditions for genes subseting, measurements (M) in rows;
            conditions table dimensions are M*D (D are description types).
            Rows should have same order  as genes table  columns (specified upon initialisation, dimension G(genes)*M) -
            column names of gene table and specified column in conditions table should match.
        :param conditions_names_column: conditions table column that matches genes index - tests that have same order
        """
        NeighbourCalculator.check_numeric(genes)
        if remove_zero:
            genes = genes[(genes != 0).any(axis=1)]
        self._genes = genes
        if conditions is not None:
            if list(conditions.loc[:, conditions_names_column]) != list(self._genes.columns):
                raise ValueError('Conditions table row order must match genes table column order.')
        self.conditions = conditions

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
                   batches: list = None, remove_batch_zero: bool = True, return_neigh_dist: bool = False,
                   genes_query_names: list = None, remove_self: bool = False, metric: str = 'cosine'):
        """
        Calculates neighbours of genes on whole gene data or its subset by column.
        :param n_neighbours: Number of neighbours to obtain for each gene
        :param inverse: Calculate most similar neighbours (False) or neighbours with inverse profile (True)
        :param scale: Scale expression by gene with 'minmax' (min=0, max=1) or 'mean0std1' (mean=0, std=1) or 'none'
        :param log: Should expression data be log2 transformed
        :param batches: Should comparisons be made for each batch separately.
            Batches should be a list of batch group names for each column (eg. length of batches is n columns of genes).
        :param remove_batch_zero: Remove genes that have all expression values 0 for each batch individually.
        :param return_neigh_dist: Instead of parsed dictionary return tuple with NN matrix and similarity matrix,
        as returned by pynndescent (converted to similarities) but named with gene names in data frame.
        :param genes_query_names: Use only the specified genes as query.
        :param remove_self: Used only if return_neigh_dist is true. Whether to remove sample from its closest
        neighbours or not. If retunr_neigh_dist is False this is done automatically. This also removes last
        column/neighbours is self is not present - should not be used with inverse.
        :param metric: cosine or correlation (pearson)
        :return: Dict with gene names as tupple keys (smaller by alphabet is first tuple value) and
            values representing cosine similarity. If batches are used such dicts are returned for each batch
            in form of dict with batch names as keys and above mentioned dicts as values. Or see return_neigh_dist.
        """
        if scale not in self.SCALES:
            raise ValueError('Scale must be:', self.SCALES)

        genes = self._genes

        if batches is None:
            if genes_query_names is not None:
                genes_query = genes.loc[genes_query_names, :]
            else:
                genes_query = None
            return NeighbourCalculator.calculate_neighbours(genes=genes, n_neighbours=n_neighbours, inverse=inverse,
                                                            scale=scale,
                                                            log=log, return_neigh_dist=return_neigh_dist,
                                                            genes_query_data=genes_query, remove_self=remove_self,
                                                            metric=metric)
        else:
            batch_groups = set(batches)
            batches = np.array(batches)
            results = dict()
            for batch in batch_groups:
                genes_sub = genes.T[batches == batch].T
                if remove_batch_zero:
                    genes_sub = genes_sub[(genes_sub != 0).any(axis=1)]

                if genes_query_names is not None:
                    genes_query_sub = genes_sub.loc[genes_query_names, :]
                else:
                    genes_query_sub = None

                result = NeighbourCalculator.calculate_neighbours(genes=genes_sub, n_neighbours=n_neighbours,
                                                                  inverse=inverse,
                                                                  scale=scale, log=log, description=batch,
                                                                  return_neigh_dist=return_neigh_dist,
                                                                  genes_query_data=genes_query_sub,
                                                                  remove_self=remove_self,
                                                                  metric=metric)
                results[batch] = result
            return results

    @staticmethod
    def calculate_neighbours(genes, n_neighbours: int, inverse: bool, scale: str, log: bool,
                             description: str = '', return_neigh_dist: bool = False,
                             genes_query_data: pd.DataFrame = None, remove_self: bool = False, metric: str = 'cosine'):
        """
        Calculate neighbours of genes.
        :param genes: Data frame as in init, gene names (rows) should match the one in init
        :param n_neighbours: Number of neighbours to obtain for each gene
        :param inverse: Calculate most similar neighbours (False) or neighbours with inverse profile (True)
        :param scale: Scale expression by gene with 'minmax' (min=0, max=1) or 'mean0std1' (mean=0, std=1) or 'none'
        :param log: Should expression data be log2 transformed
        :param description: If an error occurs in KNN index formation report this with error
        :param return_neigh_dist: Instead of parsed dictionary return tuple with NN matrix and similarities matrix,
        as returned by pynndescent (converted to similarities) but named with gene names  in data frame.
        :param genes_query_data: Use this as query. If None use genes.
        :param remove_self: Used only if return_neigh_dist is true. Whether to remove sample from its closest
        neighbours or not. If retunr_neigh_dist is False this is done automatically. This also removes last
        column/neighbours is self is not present - should not be used with inverse.
        :param metric: 'cosine' or Pearson 'correlation'
        :return: Dict with gene names as tuple keys (smaller by alphabet is first tuple value) and
            values representing cosine similarity. Or see return_neigh_dist
        """
        genes_index, genes_query = NeighbourCalculator.get_index_query(genes=genes, inverse=inverse, scale=scale,
                                                                       log=log,
                                                                       genes_query_data=genes_query_data)
        # Can set speed-quality trade-off, default is ok
        try:
            index = NNDescent(genes_index, metric=metric, n_jobs=4)
        except ValueError:
            try:
                index = NNDescent(genes_index, metric=metric, tree_init=False, n_jobs=4)
                warnings.warn(
                    'Dataset ' + description + ' index computed without tree initialisation',
                    Warning)
            except ValueError:
                raise ValueError('Dataset ' + description + ' can not be processed by pydescent')
        neighbours, distances = index.query(genes_query.tolist(), k=n_neighbours)

        if genes_query_data is None:
            genes_query_data = genes
        if return_neigh_dist:
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
        Get gene data scaled to be index or query for neighbour search.
        :param genes: Gene data for index and query.
        :param inverse: Inverse query to compute neighbours with opposite profile. True if use inverse.
        :param scale: Scale expression by gene with 'minmax' (min=0, max=1) or 'mean0std1' (mean=0, std=1) or 'none'.
        :param log: Should expression data be log2 transformed.
        :param genes_query_data: Genes data for query, if None uses genes
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
        Transform lists of neighbours and distances into dictionary with neighbours as keys and values as similarities.
        If pair of neighbours is given more than once it is overwritten the second time it is added to dictionary.
        For cosine similarity the above should be always the same.
        The neighbours (in form of index) are named based on gene data (index) stored in NeighbourCalculator instance.
        :param neighbours: Array of shape: genes*neighbours, where for each gene there are specified neighbours
        :param distances: Array of shape: genes*neighbours, where for each gene there are specified distances
            for each neighbour, as they are given above
        :param genes_query: DataFrame used to name genes in neighbours & distanmces rows (querry),
         gene names in rows.
        :param genes_idx: DataFrame used to name genes in neighbours & distanmces table entries,
         gene names in rows.
        :return: Dict with gene names as tuple keys (smaller by alphabet is first tuple value) and
            values representing similarity: 1-distance
        """
        parsed = dict()
        for gene in range(distances.shape[0]):
            for neighbour in range(distances.shape[1]):
                distance = distances[gene, neighbour]
                # Because of rounding the similarity may be slightly above one and distance slightly below 0
                if distance < 0 or distance > 1:
                    if round(distance, 4) != 0 or distance > 1:
                        warnings.warn(
                            'Odd cosine distance at ' + str(gene) + ' ' + str(neighbour) + ' :' + str(distance),
                            Warning)
                    distance = 0
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
                        # Can do average directly as there will not be more than 2 pairs with same elements
                        # (eg. both possible positions: a,b and b,a for index,query)
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
        if (np.around(distances, 4).any() < 0):
            warnings.warn(
                'Odd cosine distance in the matrix', Warning)
        return NeighbourCalculator.cosine_dist_to_sim(distances)

    @staticmethod
    def remove_self_pynn_matrix(neighbours: pd.DataFrame, similarities: pd.DataFrame):
        """
        Parse similarities and neighbours data frames. Remove entries for each gene that represent itself being
        the closest neighbour. If row name not in neghours row removes last element from neighbours and similarities.
        :param similarities: Similarities matrix from neighbours function (row names with genes)
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
        Works on cosine and Pearson correlation distances
        """
        return 1 - dist


def merge_genes_conditions(genes: pd.DataFrame, conditions: pd.DataFrame, matching) -> pd.DataFrame:
    """
    Merge dataframes with genes and conditions
    :param genes: Expression data, genes in rows, measurements in columns, dimensions G*M
    :param conditions: Description (columns) of each measurements (rows), dimensions M*D
    :param matching: Which column in conditions matches column names in genes
    :return: Data frame with merged genes and conditions
    """
    conditions = conditions.copy()
    conditions.index = conditions[matching]
    return pd.concat([genes.T, conditions], axis=1, sort=True)


def split_data(data: pd.DataFrame, split_by: str) -> dict:
    """
    Split data by column
    :param data: Data to be split by values of a column
    :param split_by: Column name for splitting
    :return: Key: split_by column value, value: data of this split_by column value
    """
    data_splitted = {}
    groupped = data.groupby(by=split_by)
    for group in groupped.groups.keys():
        data_splitted[group] = (groupped.get_group(group))
    return data_splitted


# *******************
# **** Manage data (load data, specify saving path)
# Path to expression data
path_data = '/home/karin/Documents/timeTrajectories/data/RPKUM/combined/'
path_results = '/home/karin/Documents/git/baylor_dicty_paper/'

# Load expression data
genes = pd.read_csv(path_data + 'mergedGenes_RPKUM.tsv', sep='\t', index_col=0)
conditions = pd.read_csv(path_data + 'conditions_mergedGenes.tsv', sep='\t', index_col=None)

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
    threshold_data[strain]= threshold

save_pickle(path_results + 'strainThresholds_k2_m0s1log_highest' + str(1 - THRESHOLD_QUANTILE) + '.pkl',threshold_data)

# ***********************
# *** Find close gene neighbours in each strain

# Find samples of each strain
batches = list(conditions['Strain'])
batches = np.array(batches)

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

# Load by strain similarity thresholds
threshold_data = load_pickle(path_results + 'strainThresholds_k2_m0s1log_highest' + str(1 - THRESHOLD_QUANTILE) +
                             '.pkl')

# Extract genes from strains and merge into a single matrix
files = [f for f in glob.glob(path_results + 'neighbours_kN' + str(KNN) + '_mean0std1_log_' + "*.pkl")]
n_genes = genes.shape[0]
genes_dict = dict(zip(genes.index, range(n_genes)))
merged_results = np.zeros((n_genes, n_genes))
for f in files:
    strain = f.split('_')[-1].replace('.pkl', '')
    result = load_pickle(f)
    similarity_threshold = threshold_data[strain]
    print('Adding data of',strain, ' using similarity threshold', similarity_threshold)
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
ration_max = 0.99
proportion = 0.1
threshold = genes.quantile(q=ration_max, axis=1) * proportion

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
                               '_minExpressed' + str(ration_max) + str(proportion) +
                               '_strainFilterMin' + str(min_strains) + 'Max' + str(max_strains) + '.tsv', sep='\t')

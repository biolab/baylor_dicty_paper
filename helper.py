import pandas as pd
import pickle


GROUPS = {'amiB': 'agg-', 'mybB': 'agg-', 'acaA': 'agg-', 'gtaC': 'agg-',
          'gbfA': 'lag_dis', 'tgrC1': 'lag_dis', 'tgrB1': 'tag_dis', 'tgrB1C1': 'tag_dis',
          'tagB': 'tag', 'comH': 'tag',
          'ecmARm': 'cud', 'gtaI': 'cud', 'cudA': 'cud', 'dgcA': 'cud', 'gtaG': 'cud',
          'AX4': 'WT', 'MybBGFP': 'WT',
          'acaAPkaCoe': 'sFB', 'ac3PkaCoe': 'sFB',
          'pkaR': 'prec', 'PkaCoe': 'prec'}
STAGES = ['no_agg', 'stream', 'lag', 'tag', 'tip', 'slug', 'mhat', 'cul', 'FB','yem']

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

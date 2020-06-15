print('Preparing PC1 vs time plots.')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import preprocessing as pp
from sklearn.decomposition import PCA
from sklearn.model_selection import LeaveOneOut
import itertools
import matplotlib
from collections import defaultdict
import random
import matplotlib.patches as mpatches
from matplotlib import rcParams
import os

import pygam

from helper import save_pickle, GROUPS, STAGES, PATH_RESULTS, PATH_DATA

# ***********************
# **** Helper functions


class CustomScaler:

    def __init__(self, reference: np.ndarray):
        """
        :param reference: np.ndarray or pandas DF on which to fit scaler
        """
        if isinstance(reference, pd.DataFrame):
            reference = reference.values
        self.reference = reference
        self.scalers = {}

    def transform(self, data: np.ndarray, log, scale: str):
        """

        :param data: Data to be scaled. np.ndarray od pandas DF.
        :param log: log2(data+1)
        :param scale: 'minmax','m0s1','divide_mean'
        :return: Scaled data
        """
        if not (log, scale) in self.scalers.keys():
            scaler = None
            ref = self.reference.copy()
            if log:
                ref = np.log2(ref + 1)
            if scale == 'minmax':
                scaler = pp.MinMaxScaler()
                scaler.fit(ref)
            elif scale == 'm0s1':
                scaler = pp.StandardScaler()
                scaler.fit(ref)
            elif scale == 'divide_mean':
                scaler = ref.mean(axis=0)
            self.scalers[(log, scale)] = scaler

        scaler = self.scalers[(log, scale)]
        scaled = None
        if isinstance(data, pd.DataFrame):
            data = data.values

        if log:
            data = np.log2(data + 1)
        if scale in ['minmax', 'm0s1']:
            scaled = scaler.transform(data)
        elif scale == 'divide_mean':
            scaled = data / scaler
        return scaled


def get_dimredplot_param(data, col, default, to_mode=False):
    """
    Based on data DF (information for single/multiple points) find series or mode of the parameter or use default
    if there is no column for the parameter.
    :param data: DF with data
    :param col: Column for which to extract the value
    :param default: Default to use if column is absent for the data.
    :param to_mode: Convert result to mode insetad of returning extracted series.
    :return: Extracted data as column or mode.
    """
    if isinstance(data, pd.DataFrame):
        if col in data.columns:
            result = data[col]
            if to_mode:
                result = result.mode()[0]
            return result
        else:
            return default
    elif isinstance(data, pd.Series):
        if col in data.index:
            return data[col]
        else:
            return default
    else:
        return default


# Jitter function
def rand_jitter(n, min, max, strength=0.005):
    """
    Number is jittered based on: n + randomN * (max-min) * stength, where -1 <= randomN <= 1
    :param n: Number to jitter
    :param min: Used to determine the size of the jittering
    :param max: Used to determine the size of the jittering
    :param strength: Larger strength leads to stronger jitter. Makes sense to be below 1 (adds random number scaled by
    max-min and strength.
    :return: New number.
    """
    dev = (max - min) * strength
    return n + random.uniform(-1, 1) * dev


# Colours for plotting
COLOURS_GROUP = {'agg-': '#d40808', 'lag_dis': '#e68209', 'tag_dis': '#ffb13d', 'tag': '#d1b30a', 'cud': '#4eb314',
                 'WT': '#0fa3ab', 'sFB': '#525252', 'prec': '#7010b0'}
COLOURS_STAGE = {'NA': '#d9d9d9', 'no_agg': '#ed1c24', 'stream': '#985006',
                 'lag': '#f97402', 'tag': '#d9d800', 'tip': '#66cf00', 'slug': '#008629', 'mhat': '#00c58f',
                 'cul': '#0ff2ff', 'FB': '#00b2ff', 'yem': '#666666'}


def dim_reduction_plot(plot_data: pd.DataFrame(), plot_by: str, fig_ax: tuple, order_column, colour_by_phenotype=False,
                       add_name=True, colours: dict = COLOURS_GROUP, colours_stage: dict = COLOURS_STAGE,
                       legend_groups='lower left', legend_phenotypes='upper right', fontsize=6,
                       plot_order: list = None, plot_points: bool = True, add_avg: bool = False, add_sem: bool = False,
                       sem_alpha: float = 0.1, alternative_lines: dict = None, sep_text: tuple = (30, 30),
                       phenotypes_list: list = STAGES, plot_lines: bool = True, jitter_all: bool = False,
                       point_alpha=0.5, point_size=5, jitter_strength: tuple = (0.005, 0.005)):
    """
    Plots PC1 vs time of strains, phenotype groups, and developmental stages.
    For plotting parameters that are not for individual points (e.g. line width, alpha)
    uses mode when plotting lines and legend.
    :param plot_data: Data of individual 'points'. Must have columns: 'x','y', 'Group' (for colouring),
        order_column (how to order points in line),
        and a column matching the plot_by parameter (for line plotting and names).
        Can have 'size' (point size - default param point_size),  'width' (line width),
        'alpha' (for plotting - default for points is point_alpha), 'linestyle', 'shape' (point shape),
        and phenotypes columns (matching phenotypes_list, valued 0 (absent) or 1 (present)).
    :param plot_by: Plot lines and text annotation based on this column.
    :param fig_ax: Tuple  (fig,ax) with plt elements used for plotting
    :param order_column: Order of plotting of groups from plot_by, first is plotted first.
    :param colour_by_phenotype: Whether to colours samples by phenotype. If false colour by 'Group' colour.
    :param add_name: Add text with name from plot_by groups.
    :param colours: Colours for 'Group'. Key: 'Group' value, value: colour.
    :param colours_stage: Colours for plotting stages, used if colour_by_phenotype=True. Dict with keys: from
        phenotypes_list and value: colour.
    :param legend_groups: Position for Group legend, if None do not plot
    :param legend_phenotypes: Position for stages legend, if None do not plot
    :param  fontsize: Fontsize for annotation.
    :param plot_order: Plot lines and SEMs in this order. Matching groups from plot_by.
    :param plot_points: Whether to plot points.
    :param add_avg: Average points with same x value (used for plotting lines and text positioning).
    :param add_sem: Plot SEM zones.
    :param sem_alpha: Alpha for SEM zone.
    :param phenotypes_list: List of phenotypes used to find stage columns in plot_data
    :param alternative_lines: Plot different lines than based on data points from plot_data.
        Dict with keys being groups obtained by plot_by and values tuple of lists: ([xs],[ys]). Use this also for
        text annotation.
    :param plot_lines: Whether to plot lines. Lines are plotted between points (rows) in plot_data, ordered by
        order_column
    :param jitter_all: Jitter all points. Else jitter only when multiple stages are annotated to same point, not
        jittering the first stage.
    :param point_alpha: Default alpha for points used if alpha column is absent
    :param point_size: Default size for points used if size column is absent
    :param jitter_strength: Tuple (strength_x, strength_y) used to jitter points. Use floats << 1 - based on data range.
        Higher uses more jittering.
    :param sep_text: Separate text annotations so that they do not overlap. Smaller number increases the separation
        of text annotations.  Tuple with (x,y), where x,y denote values for x and y axis.
    """
    # Sort data in order to be plotted
    if plot_order is not None:
        plot_data = plot_data.loc[
            plot_data[plot_by].map(dict(zip(plot_order, range(len(plot_order))))).sort_values().index]
    else:
        plot_order = plot_data[plot_by].unique()

    fig, ax = fig_ax

    # Plot data points
    if plot_points:
        # Either add one point per measurment (coloured by group) or multiple jitter points coloured by phenotypes
        if not colour_by_phenotype:
            for row_name, point in plot_data.iterrows():
                ax.scatter(point['x'], point['y'], s=get_dimredplot_param(point, 'size', point_size),
                           c=colours[point['Group']], alpha=get_dimredplot_param(point, 'alpha', point_alpha, True),
                           marker=get_dimredplot_param(point, 'shape', 'o', True))
        # By phenotypes
        else:
            min_x = plot_data['x'].min()
            min_y = plot_data['y'].min()
            max_x = plot_data['x'].max()
            max_y = plot_data['x'].max()
            for point in plot_data.iterrows():
                point = point[1]
                phenotypes = point[phenotypes_list]

                x = point['x']
                y = point['y']
                if jitter_all:
                    x = rand_jitter(n=x, min=min_x, max=max_x, strength=jitter_strength[0])
                    y = rand_jitter(n=y, min=min_y, max=max_y, strength=jitter_strength[1])
                # jitter when needed - e.g. multiple stages are annotated to a sample
                if phenotypes.sum() < 1:
                    ax.scatter(x, y, s=get_dimredplot_param(point, 'size', point_size),
                               c=colours_stage['NA'],
                               alpha=get_dimredplot_param(point, 'alpha', point_alpha, True),
                               marker=get_dimredplot_param(point, 'shape', 'o', True))
                elif phenotypes.sum() == 1:
                    phenotype = phenotypes[phenotypes > 0].index[0]
                    ax.scatter(x, y, s=get_dimredplot_param(point, 'size', point_size),
                               c=colours_stage[phenotype],
                               alpha=get_dimredplot_param(point, 'alpha', point_alpha, True),
                               marker=get_dimredplot_param(point, 'shape', 'o', True))
                else:
                    # Do not jitter the stage (point) of a sample that will be plotted first
                    first = True
                    for phenotype in phenotypes_list:
                        if phenotypes[phenotype] == 1:
                            if not first:
                                x = rand_jitter(n=x, min=min_x, max=max_x, strength=jitter_strength[0])
                                y = rand_jitter(n=y, min=min_y, max=max_y, strength=jitter_strength[1])
                            ax.scatter(x, y, s=get_dimredplot_param(point, 'size', point_size),
                                       c=colours_stage[phenotype],
                                       alpha=get_dimredplot_param(point, 'alpha', point_alpha, True),
                                       marker=get_dimredplot_param(point, 'shape', 'o', True))
                            first = False

    # Group for lines/names
    grouped = plot_data.groupby(plot_by)

    # Add SEM zones
    if add_sem:
        for name in plot_order:
            data_rep = grouped.get_group(name).sort_values(order_column)
            group = data_rep['Group'].values[0]
            grouped_x = data_rep.groupby(['x'])
            x_line = grouped_x.mean().index
            y_line = grouped_x.mean()['y']
            sem = grouped_x.sem()['y']
            ax.fill_between(x_line, y_line - sem, y_line + sem, alpha=sem_alpha, color=colours[group])

    # Add line between replicates' measurments - either lines connecting points, averages between points,
    # or predefined lines
    if plot_lines:
        for name in plot_order:
            data_rep = grouped.get_group(name).sort_values(order_column)
            group = data_rep['Group'].values[0]

            if alternative_lines is None:
                if not add_avg:
                    x_line = data_rep['x']
                    y_line = data_rep['y']
                else:
                    grouped_x = data_rep.groupby(['x'])
                    x_line = grouped_x.mean().index
                    y_line = grouped_x.mean()['y']
            else:
                x_line, y_line = alternative_lines[data_rep[plot_by].values[0]]

            ax.plot(x_line, y_line, color=colours[group],
                    alpha=get_dimredplot_param(data_rep, 'alpha', 0.5, True),
                    linewidth=get_dimredplot_param(data_rep, 'width', 0.5, True),
                    linestyle=get_dimredplot_param(data_rep, 'linestyle', 'solid', True))

    # Add replicate name
    if add_name:
        used_text_positions = pd.DataFrame(columns=['x', 'y'])
        x_span = plot_data['x'].max() - plot_data['x'].min()
        y_span = plot_data['y'].max() - plot_data['y'].min()
        for name in plot_order:
            data_rep = grouped.get_group(name).sort_values(order_column)
            group = data_rep['Group'].values[0]
            idx = -1
            # Add name near the line
            if alternative_lines is None:
                if not add_avg:
                    x_values = data_rep['x'].values
                    y_values = data_rep['y'].values
                else:
                    grouped_x = data_rep.groupby(['x'])
                    x_values = grouped_x.mean().index.values
                    y_values = grouped_x.mean()['y'].values
            else:
                x_values, y_values = alternative_lines[data_rep[plot_by].values[0]]
            # Make sure that names are separated enough
            x = float(x_values[idx]) + x_span / 500
            y = float(y_values[idx]) + y_span / 500
            while ((abs(used_text_positions['x'] - x) < (x_span / sep_text[0])).values &
                   (abs(used_text_positions['y'] - y) < (y_span / sep_text[1])).values).any():
                idx -= 1
                x = float(x_values[idx]) + x_span / 500
                y = float(y_values[idx]) + y_span / 500
            used_text_positions = used_text_positions.append({'x': x, 'y': y}, ignore_index=True)
            ax.text(x, y, data_rep[plot_by].values[0], fontsize=fontsize, color=colours[group])

    # Legends for groups and phenotypes
    alpha_legend = get_dimredplot_param(plot_data, 'alpha', 0.5)
    if isinstance(alpha_legend, pd.Series):
        alpha_legend = alpha_legend.median()
    if legend_groups is not None:
        patchList = []
        for name, colour in colours.items():
            data_key = mpatches.Patch(color=colour, label=name, alpha=alpha_legend)
            patchList.append(data_key)
        title = 'Phenotype'
        if colour_by_phenotype:
            title = title + ' (line)'
        legend_groups = ax.legend(handles=patchList, title=title, loc=legend_groups)

    if colour_by_phenotype and legend_phenotypes is not None:
        patchList = []
        for name, colour in colours_stage.items():
            data_key = mpatches.Patch(color=colour, label=name, alpha=alpha_legend)
            patchList.append(data_key)
        legend_stages = ax.legend(handles=patchList, title="Stage (point)", loc=legend_phenotypes)

    if legend_groups is not None:
        ax.add_artist(legend_groups)


# Use this to make sure that enforcing MIN/MAX X/Y will not cut off parts of plots (e.g. due to jittering)
def adjust_axes_lim(ax: plt.Axes, min_x_thresh: float, max_x_thresh: float, min_y_thresh: float, max_y_thresh: float):
    """
    Adjust ax limit so that it will be at least as small as min threshold and as big as max threshold.
    If min threshold is larger than existing min axes value it will not change (and vice versa for max).
    Thus the axes should be beforehand adjusted not to include padding around plot elements, as this will be
    included in min/max axes value as well.
    :param ax: Adjust range on axes object
    :param min_x_thresh: ax x_min must be at least that small.
    :param max_x_thresh: ax x_max must be at least that big.
    :param min_y_thresh: ax y_min must be at least that small.
    :param max_y_thresh: ax y_max must be at least that big.
    """
    y_min, y_max = ax.get_ylim()
    x_min, x_max = ax.get_xlim()

    if round(y_min, 3) >= round(min_y_thresh, 3):
        y_min = min_y_thresh
    else:
        print('min y was set to', y_min, 'instead of', min_y_thresh)
    if round(y_max, 3) <= round(max_y_thresh, 3):
        y_max = max_y_thresh
    else:
        print('max y was set to', y_max, 'instead of', max_y_thresh)

    if round(x_min, 3) >= round(min_x_thresh, 3):
        x_min = min_x_thresh
    else:
        print('min x was set to', x_min, 'instead of', min_x_thresh)
    if round(x_max, 3) <= round(max_x_thresh, 3):
        x_max = max_x_thresh
    else:
        print('max x was set to', x_max, 'instead of', max_x_thresh)

    ax.set_ylim(y_min, y_max)
    ax.set_xlim(x_min, x_max)


# *****************
# *** Load data
path_save = PATH_RESULTS + 'PC1vsTime/'
if not os.path.exists(path_save):
    os.makedirs(path_save)

genes = pd.read_csv(PATH_DATA + 'mergedGenes_RPKUM.tsv', sep='\t', index_col=0)
conditions = pd.read_csv(PATH_DATA + 'conditions_mergedGenes.tsv', sep='\t', index_col=None)

font = 'Arial'
matplotlib.rcParams.update({'font.family': font})

# In each strain group use different linetype for each strain
linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0, (5, 5))]
strain_linestyles = dict()
used_linestyles = defaultdict(set)
for strain in conditions['Strain'].unique():
    group = GROUPS[strain]
    for style in linestyles:
        if style not in used_linestyles[group]:
            used_linestyles[group].add(style)
            strain_linestyles[strain] = style
            break

# Default plot parameters
linewidth_mutant = 2
alpha_mutant = 0.7
linewidth_AX4 = 5
alpha_AX4 = 1
size_mutant = 30
scale_size_AX4 = linewidth_AX4 / linewidth_mutant

# *** PCA fit

# Data pre-processing parameters
LOG = True
SCALE = 'm0s1'

# Reference data (AX4 RPKUM data with non-zero expressed genes)
genes_data = genes[conditions.query('Strain =="AX4"')['Measurment']].copy()
genes_data = genes_data[(genes_data != 0).any(axis=1)].T
DATA_REFERENCE, SCALER = (genes_data, CustomScaler(genes_data))

# Scale reference
DATA_REFERENCE = pd.DataFrame(SCALER.transform(DATA_REFERENCE, log=LOG, scale=SCALE),
                              index=DATA_REFERENCE.index, columns=DATA_REFERENCE.columns)
# PCA on reference data
pca = PCA(n_components=1, random_state=0)
pca = pca.fit(DATA_REFERENCE)
save_pickle(path_save + 'PCA_AX4NonNullGenes_scale' + SCALE + 'log' + str(LOG) + '.pkl', pca)

# Use AX4-trained PCA to transform data of other strains
data_strains = genes[conditions['Measurment']].T[DATA_REFERENCE.columns]
data_strains = pd.DataFrame(SCALER.transform(data_strains, log=LOG, scale=SCALE), index=data_strains.index,
                            columns=data_strains.columns)
PCA_TRANSFORMED = pca.transform(data_strains).ravel()

# All strains PCA data for plotting
DATA_TRANSFORMED = pd.DataFrame({'y': PCA_TRANSFORMED,
                                 'linestyle': [strain_linestyles[strain] for strain in conditions['Strain']],
                                 'width': [linewidth_mutant if strain != 'AX4' else linewidth_AX4 for strain in
                                           conditions['Strain']],
                                 'alpha': [alpha_mutant if strain != 'AX4' else alpha_AX4 for strain in
                                           conditions['Strain']]
                                 })
DATA_TRANSFORMED[['x', 'Group', 'Strain', 'Replicate'] + STAGES] = conditions[
    ['Time', 'Group', 'Strain', 'Replicate'] + STAGES]
DATA_TRANSFORMED = DATA_TRANSFORMED.sort_values('x')

# *** GAM fitting to PC1 (Y) vs time (X) data

# CV parameter combinations
param_combinations = list(itertools.product(*[list(np.logspace(-6, 0, 11, base=2)), [10, 15, 20]]))
param_combinations = [{'lam': lam, 'n_splines': n_splines} for lam, n_splines in param_combinations]

# Select best GAM  for each strain and get data for plotting
strain_GAMs = dict()
for strain in conditions['Strain'].unique():
    data_transformed = DATA_TRANSFORMED.query('Strain =="' + strain + '"')
    data_transformed = data_transformed.sort_values('x')

    # CV to select GAM parameters (regularisation, n splines)
    # CV with loo x (time); e.g. all points at that x are used for testing
    splitter = LeaveOneOut()
    squarred_errors = defaultdict(list)
    x_list = data_transformed['x'].unique()
    for train_index, test_index in splitter.split(x_list):
        x_train, x_test = x_list[train_index], x_list[test_index]
        data_train = data_transformed[data_transformed['x'].isin(x_train)]
        data_test = data_transformed[data_transformed['x'].isin(x_test)]

        for param_idx, params in enumerate(param_combinations):
            gam = pygam.LinearGAM(pygam.s(0, **params))
            gam.fit(data_train['x'].values.reshape(-1, 1), data_train['y'].values.reshape(-1, 1))
            prediction = gam.predict(data_test['x'].values.reshape(-1, 1))
            # SE of all points at test location
            squared_error = (data_test['y'] - prediction) ** 2
            squarred_errors[param_idx].extend(list(squared_error.values))

    # Select params with smallest MSE
    mese = pd.DataFrame()
    for param_idx, sqes in squarred_errors.items():
        me = np.nanmean(sqes)
        mese = mese.append({'param_idx': param_idx, 'mese': me}, ignore_index=True)
    best = mese.sort_values('mese').iloc[0, :]
    params_best = param_combinations[int(best['param_idx'])]
    print(strain, 'GAM parameters: ', params_best)

    # Make the model on whole dataset for plotting
    gam = pygam.LinearGAM(pygam.s(0, **params_best))
    gam.fit(data_transformed['x'].values.reshape(-1, 1), data_transformed['y'].values.reshape(-1, 1))

    xs = np.linspace(min(x_list), max(x_list), 100)
    ys = gam.predict(xs)
    strain_GAMs[strain] = (xs, ys)

save_pickle(path_save + 'strainGAMs.pkl', strain_GAMs)

# ***********
# ** Plots

# Min/max y and x axis value for plotting - takes in account GAM and PC1 transformed data; synchronised across plots
y_values = []
for gam_xy in strain_GAMs.values():
    y_values.extend(gam_xy[1])
y_values.extend(list(DATA_TRANSFORMED['y']))
MAX_Y = max(y_values)
MIN_Y = min(y_values)
range_y = MAX_Y - MIN_Y
# Add some extra space (for line thickness/points size)
MAX_Y = MAX_Y + 0.02 * range_y
MIN_Y = MIN_Y - 0.02 * range_y

MIN_X = DATA_TRANSFORMED['x'].min()
MAX_X = DATA_TRANSFORMED['x'].max()
range_x = MAX_X - MIN_X
# Add more padding to account for jittering which is not included here
MAX_X = MAX_X + 0.05 * range_x
MIN_X = MIN_X - 0.05 * range_x

# Set little ax margins so that necessary (tight) ax range of each plot can be calculated
plt.rcParams['axes.xmargin'] = 0.01
plt.rcParams['axes.ymargin'] = 0.01

# *** Plots of fitted GAMs to individual strains with shown sample points
for strain in conditions['Strain'].unique():
    data_transformed = DATA_TRANSFORMED.query('Strain =="' + strain + '"')
    data_transformed = data_transformed.sort_values('x')

    xs, ys = strain_GAMs[strain]

    fig, ax = plt.subplots()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Time')
    ax.set_ylabel('PC1')

    ax.plot(xs, ys, 'k', lw=2, alpha=0.7)
    # Colour each replicate separately
    replicates = data_transformed['Replicate']
    replicates_unique = list(replicates.unique())
    cmap = plt.get_cmap('tab10').colors[:len(replicates_unique)]
    rep_colours = dict(zip(replicates_unique, cmap))
    ax.scatter(data_transformed['x'], data_transformed['y'], c=[rep_colours[rep] for rep in replicates], alpha=0.7)
    ax.set_title(strain, fontdict={'fontsize': 13, 'fontfamily': font})
    adjust_axes_lim(ax=ax, min_x_thresh=MIN_X, max_x_thresh=MAX_X, min_y_thresh=MIN_Y, max_y_thresh=MAX_Y)

    # Put x axis tickmarks at sample times
    sampling_times = [time for time in data_transformed['x'].unique() if time % 4 == 0]
    if strain == 'gtaC':
        sampling_times = data_transformed['x'].unique()
    ax.set_xticks(sampling_times)
    ax.set_xticks(data_transformed['x'].unique(), minor=True)
    ax.tick_params(axis='x', which='minor', length=5)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(path_save + 'GAM_' + strain + '.pdf')
    plt.close()

# *** Strain plots with developmental stages annotation
# See combined stages plot for legend
matplotlib.rcParams.update({'font.size': 15})
for strain in conditions['Strain'].unique():
    fig, ax = plt.subplots(figsize=(10, 4))
    data_plot = DATA_TRANSFORMED.query('Strain =="' + strain + '"').drop('alpha', axis=1)
    data_plot['linestyle'] = ['solid'] * data_plot.shape[0]
    data_plot['width'] = [linewidth_mutant] * data_plot.shape[0]
    # Plot each replicate with different symbol shape
    replicates = list(data_plot['Replicate'].unique())
    replicates_map = dict(zip(replicates, ['o', '^', 'd', 's', 'X', '*', 'v']))
    data_plot['shape'] = [replicates_map[rep] for rep in data_plot['Replicate']]
    dim_reduction_plot(data_plot,
                       plot_by='Strain', fig_ax=(fig, ax), order_column='x',
                       colour_by_phenotype=True, legend_groups=None, legend_phenotypes=None,
                       add_name=False, fontsize=11, colours={GROUPS[strain]: 'black'},
                       add_avg=True, alternative_lines=strain_GAMs, sep_text=(15, 30), jitter_all=True,
                       point_alpha=0.5, point_size=40, plot_lines=True, jitter_strength=(0.03, 0.05))

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Time')
    ax.set_ylabel('PC1')
    adjust_axes_lim(ax=ax, min_x_thresh=MIN_X, max_x_thresh=MAX_X, min_y_thresh=MIN_Y, max_y_thresh=MAX_Y)
    fig.suptitle(strain, fontsize=15, fontfamily=font)

    # Put x axis tickmarks at sample times
    sampling_times = [time for time in data_plot['x'].unique() if time % 4 == 0]
    if strain == 'gtaC':
        sampling_times = data_plot['x'].unique()
    ax.set_xticks(sampling_times)
    ax.set_xticks(data_plot['x'].unique(), minor=True)
    ax.tick_params(axis='x', which='minor', length=5)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(path_save + 'stagesGAM_' + strain + '.pdf')
    plt.close()

# *** Plot with GAM fits of all strains

# Order in which to plot strains for best visibility
strains = list(conditions['Strain'].unique())
plot_order = conditions[['Strain', 'Group']].copy()
plot_order['Group'] = pd.Categorical(plot_order['Group'],
                                     categories=['agg-', 'prec', 'WT', 'sFB', 'tag', 'cud', 'lag_dis', 'tag_dis'],
                                     ordered=True)
plot_order = list(plot_order.sort_values('Group')['Strain'].unique())
plot_order.remove('AX4')
plot_order = plot_order + ['AX4']

matplotlib.rcParams.update({'font.size': 13})
fig, ax = plt.subplots(figsize=(10, 10))
dim_reduction_plot(DATA_TRANSFORMED, plot_by='Strain', fig_ax=(fig, ax), order_column='x',
                   colour_by_phenotype=False, legend_groups='upper left',
                   add_name=True,
                   fontsize=13, plot_order=plot_order, plot_points=False,
                   alternative_lines=strain_GAMs, sep_text=(15, 30))

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Time')
ax.set_ylabel('PC1')
adjust_axes_lim(ax=ax, min_x_thresh=MIN_X, max_x_thresh=MAX_X, min_y_thresh=MIN_Y, max_y_thresh=MAX_Y)
fig.suptitle("GAM fits to PC1 vs time of all strains",
             fontdict={'fontsize': 13, 'fontfamily': font})
# Put x axis tickmarks at sample times
sampling_times = [time for time in DATA_TRANSFORMED['x'].unique() if time % 4 == 0]
ax.set_xticks(sampling_times)
ax.set_xticks(DATA_TRANSFORMED['x'].unique(), minor=True)
ax.tick_params(axis='x', which='minor', length=5)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(path_save + 'GAM_combined.pdf')
plt.close()

# *** Plot PC1 vs time of all strains with stage annotations
matplotlib.rcParams.update({'font.size': 15})
fig, ax = plt.subplots(figsize=(10, 10))
dim_reduction_plot(DATA_TRANSFORMED.drop('alpha', axis=1), plot_by='Strain', fig_ax=(fig, ax), order_column='x',
                   colour_by_phenotype=True, legend_groups=None, legend_phenotypes='upper left',
                   add_name=False, fontsize=11, plot_order=plot_order,
                   add_avg=True, alternative_lines=strain_GAMs, sep_text=(15, 30), jitter_all=True,
                   point_alpha=0.5, point_size=30, plot_lines=False, jitter_strength=(0.01, 0.01))

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Time')
ax.set_ylabel('PC1')
adjust_axes_lim(ax=ax, min_x_thresh=MIN_X, max_x_thresh=MAX_X, min_y_thresh=MIN_Y, max_y_thresh=MAX_Y)
a = fig.suptitle('PC1 vs time of all strains with annotated stages',
                 fontdict={'fontsize': 13, 'fontfamily': font})
sampling_times = [time for time in DATA_TRANSFORMED['x'].unique() if time % 4 == 0]
ax.set_xticks(sampling_times)
ax.set_xticks(DATA_TRANSFORMED['x'].unique(), minor=True)
ax.tick_params(axis='x', which='minor', length=5)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(path_save + 'stages_combined.pdf')
plt.close()

print('Preparing averaged data.')

import os
import pandas as pd
import numpy as np

from helper import merge_genes_conditions, split_data, GROUPS, STAGES, PATH_RESULTS, PATH_DATA

# *******************
# **** Manage data (load data, specify saving path)
# Path to data
path_results = PATH_RESULTS + 'averaged/'
if not os.path.exists(path_results):
    os.makedirs(path_results)

# Load expression data
genes = pd.read_csv(PATH_DATA + 'mergedGenes_RPKUM.tsv', sep='\t', index_col=0)
conditions = pd.read_csv(PATH_DATA + 'conditions_mergedGenes.tsv', sep='\t', index_col=None)

STRAIN_ORDER = pd.read_table(PATH_DATA + 'strain_order.tsv', header=None).values.ravel()
# *****************
# *** Expression averaged across replicates of strains by timepoint

# *** Average RPKUM data of strains across timepoints
# Combine expression and metadata
merged = merge_genes_conditions(genes=genes, conditions=conditions[['Time', 'Measurment', 'Strain']],
                                matching='Measurment')
splitted = split_data(data=merged, split_by='Strain')
# Average by strain
by_strain = {}
for split, data in splitted.items():
    genes_avg_stage = data.groupby('Time').mean()
    by_strain[split] = genes_avg_stage.T
# Combine data of different strains and add metadata
strains_data = []
for strain, data in by_strain.items():
    strain_data = pd.DataFrame({'Strain': [strain] * data.shape[1]}, index=data.columns).T.append(data)
    strains_data.append(strain_data)
genes_avg = pd.concat(strains_data, axis=1).T
genes_avg['Time'] = genes_avg.index
genes_avg.index = [strain + '_' + str(time) for strain, time in
                   zip(genes_avg['Strain'], genes_avg['Time'])]

# Add strain group
genes_avg['Group'] = [GROUPS[strain] for strain in genes_avg['Strain']]

# Sort the data
genes_avg['Strain'] = pd.Categorical(genes_avg['Strain'], categories=STRAIN_ORDER,
                                     ordered=True)
genes_avg = genes_avg.sort_values(['Strain', 'Time'])

# *** Scale the averaged data
genes_avg_scaled = genes_avg.copy()
genes_avg_scaled = genes_avg_scaled.drop(['Strain', 'Time', 'Group'], axis=1)
genes_avg_scaled = genes_avg_scaled.apply(pd.to_numeric)
# Linearly scale with Xth percentile (given as ratio) of each gene
percentile = 0.99
genes_avg_percentile = genes_avg_scaled.quantile(q=percentile, axis=0)
genes_avg_scaled = genes_avg_scaled - genes_avg_percentile
genes_avg_percentile = genes_avg_percentile.replace(0, 1)
genes_avg_scaled = genes_avg_scaled / genes_avg_percentile
# Bound values at max = 0.1 (removing extremely high values) to enable better colour scaling for visualisation
max_val = 0.1
genes_avg_scaled[genes_avg_scaled > max_val] = max_val
genes_avg_scaled[['Time', 'Strain', 'Group']] = genes_avg[['Time', 'Strain', 'Group']]

# Save the data
genes_avg_scaled.to_csv(path_results + 'genes_averaged_scaled_percentile' + str(percentile)[2:] +
                        '_max' + str(max_val) + '.tsv',
                        sep='\t')

# ******************
# *** Expression peak time in AX4

pattern_data = []
genes_avg_ax4 = genes_avg.query('Strain=="AX4"')
genes_avg_ax4.index = genes_avg_ax4['Time']
genes_avg_ax4 = genes_avg_ax4.T.drop(['Strain', 'Group', 'Time'])
for gene, data in genes_avg_ax4.iterrows():
    peak = data.sort_values().index[-1]
    pattern_data.append({'Gene': gene, 'Peak': peak})
pattern_data = pd.DataFrame(pattern_data)
pattern_data.to_csv(path_results + 'gene_peaks_AX4.tsv', sep='\t', index=False)


# ****************************
# **** Developmental stages averaged across replicates of a strain by timepoint

def set_infered(averaged: dict, infered_stages: list):
    """
    Change the averaged stages data by setting only some of the stages to no_image fill
    :param averaged: Stages data to be modified, stages as keys and 'yes'/'no'/'no_image' as values
    :param infered_stages: Stages to be set to 'no_image', others are set to 'no'
    """
    for col in STAGES:
        fill = 'no'
        if col in infered_stages:
            fill = 'no image'
        averaged[col] = fill


def avgsample_name(group: tuple) -> str:
    """
    Make a sample name from strain,time tuple obtained from pandas groupby
    :param group: Tuple with strain,time
    :return: strain_time
    """
    return group[0] + '_' + str(group[1])


# TOOD rename this file
conditions_originalstage = pd.read_csv(PATH_DATA + 'conditions_noninfered_mergedGenes.tsv', sep='\t',
                                       index_col=None).sort_values(['Strain', 'Time'])
# Make averaged stages data, averaging samples by strain and time
averaged_stages = pd.DataFrame(columns=STAGES)
grouped = conditions_originalstage.groupby(['Strain', 'Time'])
for group, data in grouped:
    name = avgsample_name(group)
    # Set all present stages (in the samples of the strain and timepoint) to yes and others to no.
    # If no stage is present set all the stages to no_image.
    averaged = {}
    for col in STAGES:
        pheno_data = data[col]
        if (pheno_data == 1).any():
            averaged[col] = 'yes'
        elif (pheno_data == -1).all():
            averaged[col] = 'no image'
        else:
            averaged[col] = 'no'

    # If no image is present for an averaged sample try to infer the stage/set some stages to 'no'.
    time = group[1]
    strain = group[0]
    # ac3PkaCoe has no images, but it is assumed that t=0 is no_agg - add 'no' to other stages.
    if strain == 'ac3PkaCoe':
        if time == 0:
            set_infered(averaged=averaged, infered_stages=['no_agg'])
    # gtaC has no images but is assumed that all samples are no_agg - add 'no' to other stages.
    elif strain == 'gtaC':
        set_infered(averaged=averaged, infered_stages=['no_agg'])
    else:
        if (np.array(list(averaged.values())) == 'no image').all():
            # It is assumed that all samples at t=0 are no_agg - add 'no' to other stages.
            if time == 0:
                set_infered(averaged=averaged, infered_stages=['no_agg'])
            # Set all stages that were 'no' in the previous timepoint to 'no'.
            else:
                strain_times = conditions_originalstage.query('Strain =="' + strain + '"')['Time'].unique()
                strain_times.sort()
                previous_time = int(strain_times[np.argwhere(strain_times == time) - 1])
                previous_data = averaged_stages.loc[avgsample_name((strain, previous_time)), :]
                previous_stages = list(previous_data[previous_data != 'no'].index)
                set_infered(averaged=averaged, infered_stages=previous_stages)

    averaged_stages = averaged_stages.append(pd.DataFrame(averaged, index=[name]), sort=True)
averaged_stages.index.name = 'Name'
# Order the columns based on stage order in development
averaged_stages = averaged_stages.reindex(STAGES, axis=1)
averaged_stages.to_csv(path_results + 'averageStages.tsv', sep='\t')

# *******************
# **** Expression averaged across main stages

# Metadata for samples with annotated main stage
conditions_main = conditions.query('~main_stage.isna()', engine='python')

# Average expression
genes_stage = genes[conditions_main.Measurment].copy().T
genes_stage = genes_stage.join(pd.DataFrame(conditions_main[['Strain', 'main_stage']].values,
                                            index=conditions_main['Measurment'], columns=['Strain', 'main_stage']))
genes_avg_stage = genes_stage.groupby(['Strain', 'main_stage']).mean().reset_index()

# Sort by strain and stage
genes_avg_stage['main_stage'] = pd.Categorical(genes_avg_stage['main_stage'], categories=STAGES, ordered=True)
genes_avg_stage['Strain'] = pd.Categorical(genes_avg_stage['Strain'], categories=STRAIN_ORDER, ordered=True)
genes_avg_stage = genes_avg_stage.sort_values(['Strain', 'main_stage'])
# Add metadata
genes_avg_stage['Group'] = [GROUPS[strain] for strain in genes_avg_stage['Strain']]
genes_avg_stage.index = [strain + '_' + main
                         for strain, main in genes_avg_stage[['Strain', 'main_stage']].values]

# *** Scale
genes_avg_scaled_stage = genes_avg_stage.copy()
genes_avg_scaled_stage = genes_avg_scaled_stage.drop(['Strain', 'main_stage', 'Group'], axis=1)
# Scale with Xth percentile (given as ratio) of each gene
percentile = 0.99
genes_avg_percentile_stage = genes_avg_scaled_stage.quantile(q=percentile, axis=0)
genes_avg_scaled_stage = genes_avg_scaled_stage - genes_avg_percentile_stage
genes_avg_percentile_stage = genes_avg_percentile_stage.replace(0, 1)
genes_avg_scaled_stage = genes_avg_scaled_stage / genes_avg_percentile_stage
# Set upper limit to scaled expression values for better representation on colour scale
max_val = 0.1
genes_avg_scaled_stage[genes_avg_scaled_stage > max_val] = max_val
genes_avg_scaled_stage[['Strain', 'main_stage', 'Group']] = genes_avg_stage[['Strain', 'main_stage', 'Group']]

genes_avg_scaled_stage.to_csv(path_results + 'genes_averaged_mainStage_scale_percentile' + str(percentile)[2:] +
                              '_max' + str(max_val) + '.tsv', sep='\t')

# *******************
# *** Averaged data from Nichols, et al. (2020)
genes_mb = pd.read_csv(PATH_DATA + 'mediaBuffer_RPKUM.tsv', sep='\t', index_col=0)
conditions_mb = pd.read_csv(PATH_DATA + 'mediaBuffer_mergedGenes.tsv', sep='\t', index_col=None)

# Average expression
genes_group_mb = genes_mb[conditions_mb.Measurment].copy().T
genes_group_mb = genes_group_mb.join(pd.DataFrame(conditions_mb[['Group', 'Time']].values,
                                            index=conditions_mb['Measurment'], columns=['Group', 'Time']))
averaged_mb = genes_group_mb.groupby(['Group', 'Time']).mean().reset_index()
# Sort by Group and Time
averaged_mb['Group'] = pd.Categorical(averaged_mb['Group'], categories=['buff','media'], ordered=True)
averaged_mb = averaged_mb.sort_values(['Group', 'Time'])

# Scale the data
percentile = 0.99
max_val = 0.1
genes_avg_scaled_mb = averaged_mb.copy()
genes_avg_scaled_mb = genes_avg_scaled_mb.drop(['Group', 'Time'], axis=1)
genes_avg_percentile_mb = genes_avg_scaled_mb.quantile(q=percentile, axis=0)
genes_avg_scaled_mb = genes_avg_scaled_mb - genes_avg_percentile_mb
genes_avg_percentile_mb = genes_avg_percentile_mb.replace(0, 1)
genes_avg_scaled_mb = genes_avg_scaled_mb / genes_avg_percentile_mb
# Limit at max 0.1 (to remove outliers in visualisation)
genes_avg_scaled_mb[genes_avg_scaled_mb > max_val] = max_val
genes_avg_scaled_mb[['Group', 'Time']] = averaged_mb[['Group', 'Time']]
genes_avg_scaled_mb.index = [group + '_' + str(time)
                                 for group, time in genes_avg_scaled_mb[['Group', 'Time']].values]

genes_avg_scaled_mb.to_csv(path_results + 'genesMediaBuffer_averaged_scaled_percentile' + str(percentile)[2:] +
                        '_max' + str(max_val) + '.tsv', sep='\t')

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm, matplotlib.colors, matplotlib.figure
import math

from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher
from orangecontrib.bioinformatics.geneset.utils import GeneSet
import orangecontrib.bioinformatics.go as go
from orangecontrib.bioinformatics.geneset.__init__ import load_gene_sets

ORGANISM = 44689

font = 'Arial'
matplotlib.rcParams.update({'font.family': font})


class GeneSetData:
    """
    Stores gene set with enrichment statistics
    """

    def __init__(self, gene_set: GeneSet, ontology: tuple = None, in_query: int = None):
        """
        :param gene_set: Gene set for which data should be stored
        :param ontology: To which ontology the gene set belongs
        :param in_query: N of genes from set present in query
        """
        self.gene_set = gene_set
        self.ontology = ontology
        self.pval = None
        self.padj = None
        self.in_query = in_query
        self.in_reference = None


def name_genes_entrez(gene_names: list, key_entrez: bool, organism: int = ORGANISM) -> dict:
    """
    Match Entrez ID to each gene name.
    :param gene_names: Gene names (eg. from dictyBase)
    :param organism: Organism ID
    :param key_entrez: If True Entrez IDs are returned as keys and names as values, False: vice versa
    :return: Dict of gene names and matching Entres IDs for genes that have Entrez ID (genes without Entrez ID are
        not returned)
    """
    entrez_names = dict()
    matcher = GeneMatcher(organism)
    matcher.genes = gene_names
    for gene in matcher.genes:
        name = gene.input_identifier
        entrez = gene.gene_id
        if entrez is not None:
            if key_entrez:
                entrez_names[entrez] = name
            else:
                entrez_names[name] = entrez
    return entrez_names


def convert_EID(genes: iter, name_EID: dict) -> set:
    """
    Convert gene names to EID based on dict as returned by name_genes_entrez.
    :param genes: Gene names
    :param name_EID: Dict where keys are names and values are EIDs
    :return: Set of gene EIDs
    """
    return set(name_EID[gene] for gene in genes if gene in name_EID.keys())


def get_gene_sets(gene_set_names: list, organism: str = ORGANISM, go_slims: bool = False,
                  set_sizes: tuple = (-np.inf, np.inf), reference: set = None) -> dict:
    """
    Get lists of gene sets sorted by ontology.
    :param gene_set_names: Names of ontologies for which to get gene sets (as returned by list_gene_sets)
    :param organism: Organism id
    :param go_slims: If ontology type (first element from tuples retured by list_gene_sets) is GO then output
        only gene sets that are in 'goslim_generic'
    :param set_sizes: Use only gene sets with size greater or equal to the 1st element of set_sizes and smaller or
        equal to the 2nd element of set_sizes, default is -inf,inf
    :param reference: List of gene EIDs to use as a reference set. Gene set filtering based on gene set size takes in
        account only genes present in the reference. If None no reference is used.
    :return: Dict with key being onotlogies name and values being GeneSet objects
    """
    gene_set_ontology = dict()

    # Get GO slims
    if go_slims:
        anno = go.Annotations(organism)
        anno._ensure_ontology()
        anno._ontology.set_slims_subset('goslim_generic')
        slims = anno._ontology.slims_subset

    # Load and filter gene sets
    for gene_set_name in gene_set_names:
        gene_sets = load_gene_sets(gene_set_name, str(organism))
        if reference is None:
            gene_sets = [gene_set for gene_set in gene_sets if set_sizes[0] <= len(gene_set.genes) <= set_sizes[1]]
        else:
            gene_sets = [gene_set for gene_set in gene_sets if
                         set_sizes[0] <= len(gene_set.genes & reference) <= set_sizes[1]]
        if go_slims and gene_set_name[0] == 'GO':
            gene_sets = [gene_set for gene_set in gene_sets if gene_set.gs_id in slims]
        gene_set_ontology[gene_set_name] = gene_sets
    return gene_set_ontology


def group_diff_enrichment(query_names: list, group: str, name_eid: dict, all_gene_names_eid: list, gene_sets_ontology,
                          padj: float = 0.25, min_overlap: int = None,
                          use_annotated_genes: bool = False,
                          make_enrichment_bar=False,
                          max_FE_bar=None, min_FDR_bar=None, cmap_FDR_bar='viridis', lFDR_base_bar=10):
    """
    Calculate and display gene group enrichment for gene sets.
    :param query_names: Gene group gene names.
    :param group: Name of the gene group used in the title.
    :param name_eid: Dict with key: gene name, value: EID.
    :param all_gene_names_eid: Reference EIDs
    :param gene_sets_ontology: Dict with ontology names as keys and values being lists of corresponding
        gene set objects.
    :param padj: FDR threshold - display only gene sets with FDR <= padj.
    :param min_overlap: Min gene group-gene set overlap. Display only gene sets with gene group overlap >= min_overlap.
    :param use_annotated_genes: If True use only genes (from reference and gene group) that have at
        least one gene set annotation.
    :param make_enrichment_bar: Make enrichment table with fold enrichment barplots.
    :param max_FE_bar: Max fold enrichment to display on enrichment barplot axis. Can be None to autoadjust.
        See plot_enrichment_bar.
    :param min_FDR_bar: Min FDR to use for enrichment barplot colours - all FDR that <= min_FDR_bar have same colour.
         See plot_enrichment_bar.
    :param cmap_FDR_bar: Colourmap for enrichment barplot, plt cmap name or list of colours, as in plot_enrichment_bar.
    :param lFDR_base_bar: log FDR base for enrichment barplot colour.  See plot_enrichment_bar.
    :return: List of: enrichment table, enrichment bar fig,ax tuple (if make_enrichment_bar is True).
    """
    # Gene group/query name to EID
    query_EID = convert_EID(genes=query_names, name_EID=name_eid)
    print('***  ' + group + ' selected:', len(query_names), 'with EID:', len(query_EID))

    reference_gene_eids = all_gene_names_eid.copy()
    query_eids = query_EID.copy()

    # If using only annotated genes filter query and reference EID
    if use_annotated_genes:
        gene_sets_genes = set()
        for gene_set_name, gene_sets in gene_sets_ontology.items():
            for gene_set in gene_sets:
                gene_sets_genes.update(gene_set.genes)
        reference_gene_eids = set(reference_gene_eids) & gene_sets_genes
        query_eids = set(query_eids) & gene_sets_genes

        query_annotated_ratio = 'NA'
        if len(query_EID) > 0:
            query_annotated_ratio = round(len(query_eids) / len(query_EID), 2)
        print("Genes annotated with a gene set in reference %.1f%% and group %.1f%%" % (
            (len(reference_gene_eids) / len(all_gene_names_eid)) * 100, query_annotated_ratio * 100))

    query_in_enriched = set()
    result = None
    fig, ax = None, None
    fig_bar, axs_bar = None, None
    if len(query_eids) > 0:
        # Calculate enrichment
        # For p value and padj calculation uses all gene sets that have overlap >=1; from gene_set_enrichment
        enrichment = gene_set_enrichment(query_eids, reference_EID=reference_gene_eids,
                                         padj_threshold=padj, min_overlap=min_overlap,
                                         gene_sets_ontology=gene_sets_ontology)
        # Enrichment table
        if len(enrichment) > 0:
            enrichment_display = list()
            enrichment = sorted(enrichment, key=lambda data: data.padj)
            for enriched in enrichment:
                query_in_enriched.update(enriched.gene_set.genes & query_eids)
                fold_enriched = (enriched.in_query / len(query_eids)) / (
                        enriched.in_reference / len(reference_gene_eids))
                enrichment_display.append({'Gene set': enriched.gene_set.name,
                                           'Ontology': enriched.ontology[0] + ': ' + enriched.ontology[1],
                                           'FDR': "{:.2e}".format(enriched.padj), 'N in group': enriched.in_query,
                                           # 'Set size': len(enriched.gene_set.genes),
                                           'N in ref.': enriched.in_reference,
                                           'Fold enrichment': fold_enriched})
            result = pd.DataFrame(enrichment_display)

            # Enrichment barplot
            if make_enrichment_bar:
                fig_bar, axs_bar = plot_enrichment_bar(df=result, query_n=len(query_eids), used_padj=padj,
                                                       reference_n=len(reference_gene_eids),
                                                       min_FDR=min_FDR_bar, max_FE=max_FE_bar, cmap=cmap_FDR_bar,
                                                       base_lFDR=lFDR_base_bar)

    print('Enrichment at FDR: ' + str(padj) + ' and min group - gene set overlap', str(min_overlap))
    print('N group genes in displayed gene sets:', len(query_in_enriched), 'out of', len(query_eids),
          'group genes used for enrichment calculation.')

    result = [result]
    if make_enrichment_bar:
        result.append((fig_bar, axs_bar))
    return result


def gene_set_enrichment(query_EID: set, reference_EID: set,gene_sets_ontology: dict = None,
                        padj_threshold: float = None,min_overlap: int = None):
    """
    Calculate enrichment for specified gene sets. Padj is calculated on combined results for all ontologies.
    Prior to p value calculation removes gene sets that do not overlap with query.
    :param query_EID: Query Ensembl IDs
    :param reference_EID: Reference Ensembl IDs
    :param gene_sets_ontology: Dict with keys being ontology names and values being gene set lists.
    :param padj_threshold: Do not return gene sets with padj below threshold. If None no filtering is performed.
    :param min_overlap: Do not return gene sets that have ovarlap with query below min_overlap (they are still
        included in the padj calculation). If None no filtering is performed.
    :return: List of GeneSetData objects with pval, padj, in_query, and in_reference filled in based on
        enrichment calculation. If no gene sets passed the filtering or overlapped the query returns an empty list.
    """
    # Calculate enrichment
    enriched = []
    for gene_set_name, gene_sets in gene_sets_ontology.items():
        for gene_set in gene_sets:
            intersect_query = len(gene_set.genes.intersection(query_EID))
            intersect_reference = len(gene_set.genes.intersection(reference_EID))
            if intersect_query > 0:
                result = gene_set.set_enrichment(reference=reference_EID, query=query_EID)
                data = GeneSetData(gene_set=gene_set, ontology=gene_set_name)
                data.pval = result.p_value
                data.in_query = intersect_query
                data.in_reference = intersect_reference
                enriched.append(data)
    if len(enriched) > 0:
        compute_padj(enriched)

        # Filter gene sets based on enrichment significance
        if padj_threshold is not None:
            enriched = [data for data in enriched if data.padj <= padj_threshold]
        if min_overlap is not None:
            enriched = [data for data in enriched if data.in_query >= min_overlap]
    return enriched


def compute_padj(data):
    """
    Add padj (FDR Benjamini-Hochberg) values to list of GeneSetData objects, based on their p values.
    :param data: List of GeneSetData objects
    """
    # Extract pvalues
    pvals = []
    for set_data in data:
        pvals.append(set_data.pval)
    # Add padj values
    padjvals = multipletests(pvals=pvals, method='fdr_bh', is_sorted=False)[1]
    for data, padj in zip(data, padjvals):
        data.padj = padj


def plot_enrichment_bar(df: pd.DataFrame, query_n, reference_n, used_padj, min_FDR=10 ** (-10), fig_w=15, max_FE=None,
                        base_lFDR=10, cmap='viridis', automatic_size: bool = True):
    """
    Plot enrichment table with fold enrichment barplot. The table is sorted descendingly by fold enrichment.
    :param df: Enrichment table as made in group_diff_enrichment
    :param query_n: Number of query genes used for enrichment calculation.
    :param reference_n: Number of reference genes used for enrichment calculation.
    :param used_padj: Padj threshold used for enrichment filtering, used as minimal colour threshold in barplot.
    :param min_FDR: Maximal FDR to be given the maximal colour in barplot, if FDR<=min_FDR the colour is maximal.
    :param fig_w: Figure w in inches, ignored in automatic_size is True.
    :param max_FE: Max  fold enrichment used for barplot axis.
    :param base_lFDR: Log base for FDR transformation. Transforms FDRs based on -logBASE(FDR).
    :param cmap: Colourmap; plt name or list of colours. Maximal colours are used for smaller FDR (larger -log(FDR)).
    :param automatic_size: Set automatic figsize and colwidth.
    :return: fig, ax
    """
    # Create and order table for plotting
    df_plot = pd.DataFrame()
    df_plot['colour'] = [-log_base(float(padj), base_lFDR) if -log_base(float(padj), base_lFDR) <= -log_base(
        min_FDR, base_lFDR) else -log_base(min_FDR, base_lFDR)
                         for padj in df['FDR']]
    df_plot['Fold enrichment'] = df['Fold enrichment']
    df_plot['Term'] = df['Gene set']
    df_plot['Ontology'] = [ont.replace('biological_process', 'BP').replace(
        'cellular_component', 'CC').replace(
        'molecular_function', 'MF').replace(
        'Pathways', 'Path.').replace(
        'Dictybase: Phenotypes', 'DB: Pheno.').replace(
        'Custom: Baylor', 'Custom') for ont in df['Ontology'].values]
    df_plot['Group'] = ["%d (%.1f%%)" % (n, (100 * n / query_n)) for n in df['N in group']]
    df_plot['Reference'] = ["%d (%.1f%%)" % (n, (100 * n / reference_n)) for n in df['N in ref.']]
    df_plot['FDR'] = df['FDR']
    df_plot = df_plot.sort_values('Fold enrichment', ascending=False)

    # Make the figure
    return plot_table_barh(df=df_plot, bar_col='Fold enrichment', colour_col='colour', col_widths=[37, 14, 8, 9, 9],
                           fontsize=8,
                           # Figsize height based on nrows+header, add some constant needed for the tables
                           # that have less rows. Not used if automatic_size is True.
                           figsize=(fig_w, 0.33 * (df_plot.shape[0] + 1) + 0.1),
                           max_col=-log_base(min_FDR, base_lFDR), min_col=-log_base(float(used_padj), base_lFDR),
                           min_bar=1, max_bar=max_FE,
                           format_bar_axes=math.ceil, cmap=cmap, automatic_size=automatic_size)


def log_base(x: float, base):
    """
    Calculate log with specific base
    :param x: Value to log transform.
    :param base: Log base
    :return: Log transformed value
    """
    return np.log(x) / np.log(base)


def plot_table_barh(df: pd.DataFrame, bar_col, colour_col, col_widths: list, figsize: tuple, show_barcol: bool = False,
                    show_colour_col: bool = False, fontsize: int = 10,
                    min_bar: float = None, max_bar: float = None, max_col: float = None, min_col: float = None,
                    format_bar_axes=float, cmap='viridis',
                    automatic_size: bool = True):
    """
    Plot table with one of the columns represented as horizontal bars.
    :param df: Table to plot
    :param bar_col: Column to use for bar length.
    :param colour_col: Column to use for bar colouring.
    :param col_widths: Specify column withs of the table as a ratio. Ignored if automatic_size is True.
    :param figsize: Figure size. Ignored if automatic_size is True.
    :param show_barcol: Show in table the bar_col column.
    :param show_colour_col:  Show in table the colour_col column.
    :param fontsize: Image font.
    :param min_bar: Axis min for barplot.
    :param max_bar: Axis max for barplot. If None use max of values and formated max value.
    :param max_col: Maximal value for colour scale. If val>max_col use max_col colour.
    :param min_col: Minimal value for colour scale. If val<min_col use min_col colour.
    :param format_bar_axes: Bar axes labels are formated with this function.
        E.g. use int to use integers for the axis labels.
    :param cmap: Plt colourmap name or list of colours or cmap to use in linear barplot colourmap.
    :param automatic_size: Scale the table columns and figsize based on textsize. Bar plot is set to dpi=100.
    :return: fig,ax
    """
    # Colour palette range
    if min_col is None:
        min_col = df[colour_col].min()
    if max_col is None:
        max_col = df[colour_col].max()
    col_norm = matplotlib.colors.Normalize(vmin=min_col, vmax=max_col, clip=True)
    # Mappable colours
    if isinstance(cmap, str):
        cmap = matplotlib.cm.get_cmap(cmap)
    elif isinstance(cmap, list):
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", cmap)
    elif isinstance(cmap, matplotlib.colors.Colormap):
        pass
    else:
        raise ValueError('cmap must be str or list of colours or cmap')

    # Bar plot axes length scale
    if min_bar is None:
        min_bar = df[bar_col].min()
    if max_bar is None:
        max_bar = df[bar_col].max()
        max_bar = max(max_bar, format_bar_axes(max_bar))

    # Columns to show in the table
    cols = list(df.columns)
    if not show_barcol:
        cols.remove(bar_col)
    if not show_colour_col:
        cols.remove(colour_col)
    n_row = df.shape[0]

    # Create image so that each row is a subplot and table/barplots are further subplots of a row.
    border = 0.05
    width_ratios = [4, 1]
    dpi = 100
    gs_kw = dict(width_ratios=width_ratios)
    # Automatic image and col width scaling
    if automatic_size:
        ax_ws, ax_hs = calculate_table_size(df[cols].append(pd.DataFrame(cols, index=cols).T), fontsize=fontsize,
                                            space_w=15, space_h=5)
        # Heights
        # Space for headers
        ax_table_h = sum(ax_hs)
        row_h = (ax_table_h / (n_row + 1))
        header_h = ax_hs[-1] * 2
        height_ratios = [header_h] + ax_hs[:-1]
        ax_h = row_h * n_row + header_h

        # Widths
        ax_table_w = sum(ax_ws.values())
        barplot_w = 100
        ax_w = ax_table_w + barplot_w
        width_ratios = [ax_table_w, barplot_w]

        figsize = calculate_fig_size(ax_w, ax_h, borders_w=border * 2, borders_h=border * 2, dpi=dpi)
        gs_kw = dict(width_ratios=width_ratios, height_ratios=height_ratios)
    fig, axs = plt.subplots(nrows=n_row + 1, ncols=2, sharex=True, gridspec_kw=gs_kw, figsize=figsize, dpi=dpi,
                            subplotpars=matplotlib.figure.SubplotParams(
                                left=border, bottom=border, right=1 - border, top=1 - border, wspace=0, hspace=0))
    # Remove subplots edges
    for axs_row in axs:
        for ax in axs_row:
            for spine in ax.spines.values():
                spine.set_edgecolor('gray')
            ax.spines['bottom'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.yaxis.set_visible(False)
            ax.xaxis.set_visible(False)
    fig.subplots_adjust(wspace=0, hspace=0)

    # Header
    if automatic_size:
        col_widths = [ax_ws[col] for col in cols]
    the_table = axs[0][0].table(cellText=[cols], bbox=[0, 0, 1, 1], colWidths=col_widths, fontsize=fontsize, edges='',
                                cellLoc='left')
    # Make the lower border of table header
    axs[1][0].spines['top'].set_visible(True)
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(fontsize)
    for key, cell in the_table.get_celld().items():
        cell.PAD = 0

    # Plot rows
    for idx in range(n_row):
        idx_plot = idx + 1

        # Table row
        if automatic_size:
            col_widths = [ax_ws[col] for col in cols]
        the_table = axs[idx_plot][0].table(cellText=df.iloc[idx, :][cols].values.reshape(1, -1), bbox=[0, 0, 1, 1],
                                           colWidths=col_widths, fontsize=fontsize, edges='', cellLoc='left')
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(fontsize)
        for key, cell in the_table.get_celld().items():
            cell.PAD = 0

        # Barplot
        if idx == 0:
            axs[idx_plot][1].xaxis.set_visible(True)
            axs[idx_plot][1].set_xticks(
                [format_bar_axes(min_bar), format_bar_axes((min_bar + max_bar) / 2), format_bar_axes(max_bar)])
            axs[idx_plot][1].tick_params(axis='x', labelsize=fontsize)
            axs[idx_plot][1].xaxis.set_ticks_position('top')
            axs[idx_plot][1].xaxis.set_label_position('top')
            axs[idx_plot][1].set_xlabel(bar_col, fontsize=fontsize)
            axs[idx_plot][1].spines['top'].set_visible(True)
        axs[idx_plot][1].spines['left'].set_visible(True)
        axs[idx_plot][1].barh(0, df.iloc[idx, :][bar_col], 1,
                              color=cmap(col_norm(df.iloc[idx, :][colour_col])))
        axs[idx_plot][1].set_xlim([min_bar, max_bar])

    return fig, axs


def calculate_table_size(df: pd.DataFrame, fontsize: int, space_w: int, space_h: int):
    """
    Calculate table size for plt tables with no borders/spacing/padding and no index/header. Padding can be
    added to the calculation as space_w/h.
    :param df: Table
    :param fontsize: Table fontsize
    :param space_w: Horizontal space to be added to each column's min width dimensions, in dpi
    :param space_h: Vertical space to be added to each column's min height dimensions, in dpi
    :return: W,H; both in dpi. W: dict with col names as keys and widths as values. H: List of row heights
    """
    widths = {}
    heights = []
    for col in df.columns:
        max_w = 0
        for val in df[col]:
            w, h = text_dpi(str(val), fontsize=fontsize)
            if w > max_w:
                max_w = w
        widths[col] = max_w + space_w
    for row in df.index:
        max_h = 0
        for val in df.loc[row, :]:
            w, h = text_dpi(str(val), fontsize=fontsize)
            if h > max_h:
                max_h = h
        heights.append(max_h + space_h)

    return widths, heights


def calculate_fig_size(ax_w, ax_h, borders_w, borders_h, dpi=100):
    """
    Calculate figure size in inches (for figsize plt param) based on  plotting area size.
    :param ax_w: Plotting area width in dpi
    :param ax_h: Plotting area height in dpi
    :param borders_w: Ratio of width taken up by borders (non-plotting area) in dimensions of figsize.
    :param borders_h: Ratio of height taken up by borders  (non-plotting area) in dimensions of figsize.
    :param dpi: Used dpi.
    :return: fig width, fig height in inches
    """
    figw = float(ax_w / dpi) / (1 - borders_w)
    figh = float(ax_h / dpi) / (1 - borders_h)
    return figw, figh


def text_dpi(text: str, fontsize: int):
    """
    Plt text dimensions in dpi.
    :param text: Text
    :param fontsize: Text fontsize
    :return: Text width, text height in dpi.
    """
    fig, ax = plt.subplots()
    text_plot = ax.text(0, 0, text, fontsize=fontsize)
    renderer = fig.canvas.get_renderer()
    bb = text_plot.get_window_extent(renderer=renderer)
    plt.close(fig)
    return bb.width, bb.height


def plot_legend_enrichment_bar(cmap, min_FDR: float, used_padj: float, base=10):
    """
    Plot legend for enrichment barplot. Lower FDR has 'higher' colour. Transform FDR values: -logB(FDR), B is base.
    :param cmap: As in plot_linear_legend - plt cmap name or list of colours or cmap instance.
    :param min_FDR: The upper limit for FDR colour legend - all FDRs below this have same colour.
    :param used_padj: The lower limit for FDR colour legend - the threshold used for enrichment result filtering.
    :param base: Log base for FDR transform used for colour mapping.
    :return: fig, ax of the legend
    """
    label = '-log' + str(base) + '(FDR)'
    min_col = -log_base(used_padj, base)
    max_col = -log_base(min_FDR, base)
    ticks = range(math.ceil(min_col), math.floor(max_col) + 1)
    return plot_continous_legend(cmap=cmap, min_col=min_col, max_col=max_col, ticks=ticks, label=label)


def plot_continous_legend(cmap, min_col, max_col, ticks, label):
    """
    Plot a contonious colourscale legend.
    :param cmap: Plt cmap name or list of colours or cmap object
    :param min_col: Min value mapped to cmap.
    :param max_col: Max value mapped to cmap.
    :param ticks: Legend ticks in the same units as min_col and max_col
    :param label: Legend title
    :return: fig, ax of the legend
    """
    # Prepare cmap and scaling
    col_norm = matplotlib.colors.Normalize(vmin=min_col, vmax=max_col, clip=True)
    if isinstance(cmap, str):
        cmap = matplotlib.cm.get_cmap(cmap)
    elif isinstance(cmap, list):
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", cmap)
    elif isinstance(cmap, matplotlib.colors.Colormap):
        pass
    else:
        raise ValueError('cmap must be str or list of colours')

    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=col_norm)
    sm.set_array([])

    # Plot legend
    figsize = (0.8, 2.7)
    fontsize = 10
    fig, ax = plt.subplots(figsize=figsize)
    clb = fig.colorbar(sm, cax=ax, use_gridspec=True, ticks=ticks)
    for t in clb.ax.get_yticklabels():
        t.set_fontsize(fontsize)
    clb.ax.set_title(label + '\n', fontsize=fontsize, loc='left')
    margins = {"left": 0.1, "bottom": 0.02, "right": 0.5, "top": 0.9}
    fig.subplots_adjust(**margins)
    return fig, ax

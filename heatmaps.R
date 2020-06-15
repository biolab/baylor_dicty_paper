print('Preparing expression heatmaps.')

library(ComplexHeatmap)
library(circlize)
library(viridis)
library(extrafont)
library(proxy)
library(seriation)
library(dendextend)
library(purrr)

# Import fonts
font_import(prompt = FALSE)
loadfonts(device = "postscript")
loadfonts(device = "pdf")

# ************
# ** Load data shared across heatmaps

path_data<- 'Data/'
path_save <- 'Results/heatmaps/'
path_avg <- 'Results/averaged/'
if(!dir.exists(path_save)) dir.create(path_save,recursive=TRUE)

# Strain order - single column with ordered strain names
strain_order <- as.vector(read.table(paste(path_data, "strain_order.tsv", sep = ''))[, 1])

# Some plotting parameters
phenotypes_font <- 10
legend_height <- 1.5
legend_width <- 0.7
top_annotation_height <- 0.6
phenotype_annotation_height <- 3
cluster_font <- 15
fontfamily <- 'Arial'
gap_units <- 'mm'

# Colours of strain groups
group_cols <- c('agg-' = '#ed1c24', 'lag_dis' = '#f97402', 'tag_dis' = '#ffb100', 'tag' = '#d9d800', 'cud' = '#008629',
                'WT' = '#00b2ff', 'sFB' = '#1925ae', 'prec' = '#a400d4')
group_cols_text <- c('agg-' = 'black', 'lag_dis' = 'black', 'tag_dis' = 'black', 'tag' = 'black', 'cud' = '#eeeeee',
                     'WT' = 'black', 'sFB' = '#eeeeee', 'prec' = '#eeeeee')

# Colours of phenotype annotations
phenotype_cols <- c('no image' = '#d9d9d9', 'no_agg' = '#ed1c24', 'stream' = '#985006', 'lag' = '#f97402',
                    'tag' = '#d9d800', 'tip' = '#66cf00', 'slug' = '#008629', 'mhat' = '#00c58f',
                    'cul' = '#0ff2ff', 'FB' = '#00b2ff', 'yem' = '#666666')

#' Use optically ordered hc to sort genes based on AX4 expression pattern
#' @param genes List of genes to order
#' @param avg_expression DF with genes in columns and expression samples in rows. It must hae Strain column
#'    so that AX4 data can be selected.
#' @returns Ordered list of genes
optically_order_genes <- function(genes, avg_expression) {
  if (length(genes) > 1) {
    expression <- t(avg_expression[avg_expression$Strain == 'AX4', genes])
    distances <- dist(expression, method = "Euclidean")
    hc <- hclust(d = distances, method = "ward.D2")
    hc_ordered <- reorder(x = hc, dist = distances)
    genes <- as.dendrogram(hc_ordered) %>% labels
  }
  return(genes)
}

# ********************
# *** Regulons heatmap

# *** Rename regulons so that regulon with lower n umber will be plotted first.
# Needs to be done only once as the regulons file is overwritten
path_regulons<- 'Results/regulons/'
# Expression patterns file: Expression pattern in AX4.
# Genes in the first column (used as row names); a column Peak with peak time of expression in AX4 -
# used for sorting the regulons.
expression_patterns <- read.table(paste(path_avg, "gene_peaks_AX4.tsv", sep = ''),
                                  header = TRUE, row.names = 1, sep = "\t")

#' Sort clusters based on median and mean pattern time (e.g. peak time in AX4)
#' @param regulons Regulons dataframe
#' @param expression_patterns Expression pattern data frame (described above)
#' @param pattern_type Column from expression_pattern data frame (described above) used for sorting the regulons
sort_clusters <- function(regulons, expression_patterns = expression_patterns, pattern_type = 'Peak') {
  cluster_patterns_mean <- c()
  cluster_patterns_median <- c()
  clusters<- unique(regulons$Cluster)
  # Extract median and mean pattern time for each regulon
  for (cluster in clusters) {
    genes <- as.character(regulons[regulons$Cluster == cluster, 'Gene'])
    pattern_mean <- mean(expression_patterns[genes, pattern_type])
    pattern_median <- median(expression_patterns[genes, pattern_type])
    cluster_patterns_mean <- c(cluster_patterns_mean, pattern_mean)
    cluster_patterns_median <- c(cluster_patterns_median, pattern_median)
  }
  # Sort by median and then mean
  cluster_order <- data.frame('Cluster' = clusters, 'Pattern_mean' = cluster_patterns_mean, 'Pattern_median' = cluster_patterns_median)
  cluster_order <- cluster_order[order(cluster_order$Pattern_median, cluster_order$Pattern_mean),]
  return(cluster_order)
}

for(regulons_file in c("mergedGenes_minExpressed0.990.1Strains1Min1Max18_clustersAX4Louvain0.8m0s1log.tab",
                       'mergedGenes_minExpressed0.990.1Strains1Min1Max18_clustersLouvain0.4minmaxNologPCA30kN30.tab')){
  # Load regulons
  # Regulon groups tab file: First column lists genes and
  # a column named Cluster specifying cluster/regulon of each gene
  regulons <- read.table(paste(path_data, regulons_file, sep = ''), header = TRUE, sep = "\t")
  #Name the first column (should contain genes
  colnames(regulons)[1] <- 'Gene'

  # Rename the regulons in plotting order
  cluster_order <- sort_clusters(regulons = regulons, expression_patterns = expression_patterns, pattern_type = 'Peak')
  # Rename the regulons based on the new ordering and save the renamed regulons
  cluster_map<-c(paste('C',c(1:nrow(cluster_order)),sep=''))
  names(cluster_map)<-as.vector(cluster_order$Cluster)
  remap_cluster<-function(x){return(cluster_map[[x]])}
  regulons['Cluster']<-unlist(map(as.character(regulons$Cluster),remap_cluster))

  # Save renamed regulons
  write.table(regulons,paste(path_regulons,regulons_file,sep=''),row.names=FALSE, sep="\t")
}

# *** Load data and helper functions

# Averaged expression tab file: Genes in columns (already scaled), averaged strain data in rows,
# three additional comlumns: Time, Strain, and Group (meaning strain group)
avg_expression <- read.table(paste(path_avg, "genes_averaged_scaled_percentile99_max0.1.tsv", sep = ''),
                             header = TRUE, row.names = 1, sep = "\t")

# Phenotypes tab file: Short averaged sample names in rows (as in avg_expression) and columns with phenotypes.
# Phenotypes should have values: yes, no, no image
avg_phenotype <- read.table(paste(path_avg, "averageStages.tsv", sep = ''),
                            header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
# Modify avg_phenotypes data so that each stage can be coloured differently - set yes to the stage it represents
avg_phenotype[avg_phenotype == 'no'] = NA
for (col in colnames(avg_phenotype)) {
  new_col <- avg_phenotype[col]
  new_col[new_col == 'yes'] = col
  avg_phenotype[col] <- new_col
}

strain_gap <- 1
group_gap <- 2.5

# Annotation for regulons heatmap
make_annotation <- function(phenotypes_font = parent.frame()$phenotypes_font,
                            legend_height = parent.frame()$legend_height,
                            legend_width = parent.frame()$legend_width,
                            top_annotation_height = parent.frame()$top_annotation_height,
                            phenotype_annotation_height = parent.frame()$phenotype_annotation_height,
                            cluster_font = parent.frame()$cluster_font) {

  group_data <- t(avg_expression['Group'])
  rownames(group_data) <- c('Phenotypic group')
  # Strain group colours, heatmap column blocks by strain, and heatmap column blocks gaps
  group_cols_ordered <- c()
  groups_ordered <- c()
  text_cols_ordered <- c()
  gaps <- c()
  previous_group <- NULL
  for (strain in strain_order) {
    group <- as.character(avg_expression[avg_expression$Strain == strain, 'Group'][1])
    groups_ordered <- append(groups_ordered, group)
    group_cols_ordered <- append(group_cols_ordered, group_cols[group])
    text_cols_ordered <- append(text_cols_ordered, group_cols_text[group])
    # Gaps - if previous group was different add larger gap
    if (!is.null(previous_group)) {
      if (previous_group == group) {
        gaps <- append(gaps, strain_gap)
      }else {
        gaps <- append(gaps, group_gap)
      }
    }
    previous_group <- group
  }
  gaps <- unit(gaps, 'mm')

  # Time colours
  times <- unique(avg_expression$Time)
  col_time <- colorRamp2(c(min(times), max(times)), c("white", "#440154FF"))

  # Time, strain, and strain group annotation
  ht_list <- Heatmap(
  # Time annotation
  t(avg_expression['Time']), height = unit(top_annotation_height, "cm"),
  column_split = factor(avg_expression$Strain,
                        #** Ordering of the strains (column blocks) in the heatmap
                        levels = strain_order),
  column_title = NULL, column_gap = gaps,
  cluster_columns = FALSE, show_column_names = FALSE, name = '\nTime\n', col = col_time,
  heatmap_legend_param = list(at = c(min(times),
                                     as.integer(mean(c(min(times), max(times)))), max(times)),
                              grid_width = unit(legend_width, "cm"),
                              grid_height = unit(legend_height, "cm"),
                              labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                              title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)),
  row_names_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
  top_annotation = HeatmapAnnotation(
  # Phenotype group annotation
  Phenotype = anno_block(gp =
                           gpar(fill = group_cols_ordered, col = group_cols_ordered, lwd = 2, linejoin = 'mitre'),
                         labels = groups_ordered, labels_gp = gpar(col = text_cols_ordered,
                                                                   fontsize = cluster_font, fontfamily = fontfamily),
                         show_name = TRUE),
  # Strain annotation
  Strain = anno_block(gp =
                        gpar(fill = 'white', col = group_cols_ordered, lwd = 2, linejoin = 'mitre'),
                      labels = strain_order, labels_gp = gpar(col = 'black',
                                                              fontsize = cluster_font, fontfamily = fontfamily),
                      show_name = TRUE),
  annotation_name_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
  )
  )

  # Phenotype annotation
  ht_phenotype <- Heatmap(t(avg_phenotype)[, rownames(avg_expression)], height = unit(phenotype_annotation_height, "cm"),
                          cluster_columns = FALSE, cluster_rows = FALSE, show_column_names = FALSE,
                          name = '\nMorphological \nstage\n', col = phenotype_cols,
                          row_names_gp = gpar(fontsize = phenotypes_font, fontfamily = fontfamily), na_col = "white",
                          row_title = 'Morphological stage', row_title_side = 'right',
                          row_title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                          heatmap_legend_param = list(grid_width = unit(legend_width, "cm"),
                                                      grid_height = unit(legend_height, "cm"),
                                                      labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                                                      title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
                          )
  )

  # Combined annotation
  ht_list <- ht_list %v% ht_phenotype

  return(ht_list)
}

# Regulons reference (AX4) groups tab file (used for side annotation): First column lists genes and
# a column named Cluster specifying cluster/regulon of each gene
regulons2 <- read.table(paste(path_regulons,
                              "mergedGenes_minExpressed0.990.1Strains1Min1Max18_clustersAX4Louvain0.8m0s1log.tab",
                              sep = ''),header = TRUE, sep = "\t")
rownames(regulons2) <- regulons2[, 1]
regulons2 <- regulons2[, 'Cluster', drop = F]

# Sort reference (AX4 based) regulons
regulons2_temp <- data.frame(regulons2)
regulons2_temp$Gene <- row.names(regulons2_temp)
cluster_order2 <- sort_clusters(regulons = regulons2_temp, expression_patterns = expression_patterns, pattern_type = 'Peak')

# Visually order genes within reference (AX4) regulons based on AX4 averaged scaled expression data.
AX4_ordered <- c()
for (cluster in cluster_order2$Cluster) {
  genes <- row.names(regulons2)[regulons2$Cluster == cluster]
  genes <- optically_order_genes(genes = genes, avg_expression = avg_expression)
  AX4_ordered <- append(AX4_ordered, genes)
}

# reference (AX4) regulon colours: 13 distinct colours ordered by rainbow
colours_regulons2 <- c('#800000', '#e6194b', '#f58231', '#9a6324', '#ffe119', '#9dd100', '#3cb44b', '#aaffc3',
                       '#46f0f0', '#2e97ff', '#000075', '#911eb4', '#f032e6')
# refernce (AX4) regulons colour map
colours_regulons2_map <- colours_regulons2[1:length(unique(regulons2$Cluster))]
regulons2_clusters <- as.vector(unique(regulons2$Cluster))
names(colours_regulons2_map) <- regulons2_clusters[order(nchar(regulons2_clusters), regulons2_clusters)]

# *** Plot AX4 or all strains based regulons
# Plot either 'all' strains based or 'AX4' based regulons - set this for regulons_type
regulons_type <- 'all'
regulons_file <- NULL
regulon_to_letters <- NULL
if (regulons_type == 'AX4') {
  regulons_file <- "mergedGenes_minExpressed0.990.1Strains1Min1Max18_clustersAX4Louvain0.8m0s1log.tab"
  regulon_to_letters <- TRUE
} else if (regulons_type == 'all') {
  regulons_file <- 'mergedGenes_minExpressed0.990.1Strains1Min1Max18_clustersLouvain0.4minmaxNologPCA30kN30.tab'
  regulon_to_letters <- FALSE
}
# Regulon groups tab file: First column lists genes and
# a column named Cluster specifying cluster/regulon of each gene
regulons <- read.table(paste(path_regulons, regulons_file, sep = ''), header = TRUE, sep = "\t")
#Name the first column (should contain genes
colnames(regulons)[1] <- 'Gene'

# Get regulons - list unique and sort
clusters <- unique(regulons$Cluster)
vals <- as.numeric(gsub("C", "", clusters))
clusters <- clusters[order(vals)]

# Expression range for legend
expressions <- within(avg_expression, rm('Time', 'Strain', 'Group'))
min_expression <- min(expressions[, regulons$Gene])
max_expression <- max(expressions[, regulons$Gene])

# Make heatmap annotations
ht_list <- make_annotation()

# Order regulons
cluster_order <- sort_clusters(regulons = regulons, expression_patterns = expression_patterns, pattern_type = 'Peak')

# Plot concatenated regulons expression heatmaps
# Some annotations are added only for the first regulon heatmap (e.g. legend, etc.)
first <- TRUE
for (cluster in cluster_order$Cluster) {
  genes <- as.character(regulons[regulons$Cluster == cluster, 'Gene'])
  genes <- as.character(genes[order(match(genes, AX4_ordered))])
  regulons2_annotation <- rowAnnotation(AX4_clusters = regulons2[genes,], col = list(AX4_clusters = colours_regulons2_map),
                                        show_legend = FALSE, annotation_name_side = "top", show_annotation_name = first,
                                        annotation_name_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily))

  # Remove 'C' from cluster name
  # The as.character ensures that the code works with numeric clusters
  cluster_anno <- gsub('C', '', as.character(cluster))
  # Rename cluster number to a letter for AX4 regulons
  if (regulons_type == 'AX4') cluster_anno <- LETTERS[as.integer(cluster_anno)]

  heatmap <- Heatmap(t(avg_expression[, genes]), cluster_columns = FALSE, cluster_rows = FALSE, show_column_names = FALSE,
                     show_row_names = FALSE, col = viridis(256), column_title = NULL,
                     row_title = cluster_anno,
                     show_heatmap_legend = first, heatmap_legend_param = list(
      title = "\nRelative \nexpression\n",
      at = c(min_expression, round(mean(c(min_expression, max_expression)), 1), max_expression),
      grid_width = unit(legend_width, "cm"), grid_height = unit(legend_height, "cm"),
      labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily), title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)),
                     # Cluster name fontsize
                     row_title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                     left_annotation = regulons2_annotation)
  first <- FALSE
  ht_list <- ht_list %v% heatmap
}
# Plot and save the heatmap
pdf(paste0(path_save, 'expressionHeatmap_regulons_', regulons_file, '.pdf'), width = 35, height = 25)
ht_list
graphics.off()

# ***************************
# ****** Milestones heatmap

# *** Load data and helper functions
path_milestones<-'Results/milestones/'

avg_expression_stages <- read.table(
  paste0(path_avg, "genes_averaged_mainStage_scale_percentile99_max0.1.tsv"),
  header = TRUE, row.names = 1, sep = "\t")

strain_gap <- 1
group_gap <- 3.5

# Milestones heatmap annotation
make_anno_mainstage <- function(legend_height = parent.frame()$legend_height,
                                legend_width = parent.frame()$legend_width,
                                top_annotation_height = parent.frame()$top_annotation_height,
                                cluster_font = parent.frame()$cluster_font,
                                strain_gap = parent.frame()$strain_gap, group_gap = parent.frame()$group_gap,
                                gap_units = parent.frame()$gap_units) {

  group_data <- t(avg_expression_stages['Group'])
  rownames(group_data) <- c('Phenotypic group')

  # Group annotation colours and gaps
  group_cols_ordered <- c()
  groups_ordered <- c()
  text_cols_ordered <- c()
  gaps <- c()
  previous_group <- NULL
  for (strain in strain_order) {
    if (strain %in% avg_expression_stages$Strain) {
      group <- as.character(avg_expression_stages[avg_expression_stages$Strain == strain, 'Group'][1])
      groups_ordered <- append(groups_ordered, group)
      group_cols_ordered <- append(group_cols_ordered, group_cols[group])
      text_cols_ordered <- append(text_cols_ordered, group_cols_text[group])
      # Gaps - if previous group was different add larger gap
      if (!is.null(previous_group)) {
        if (previous_group == group) {
          gaps <- append(gaps, strain_gap)
        }else {
          gaps <- append(gaps, group_gap)
        }
      }
      previous_group <- group
    }else {
      strain_order <- strain_order[strain_order != strain]
    }
  }
  gaps <- unit(gaps, gap_units)

  ht_list <- Heatmap(
  # Main stage annotation
  t(avg_expression_stages['main_stage']),
  height = unit(top_annotation_height, "cm"),
  column_split = factor(avg_expression_stages$Strain,
                        # Ordering of the strain column blocks in the heatmap
                        levels = strain_order),
  column_title = NULL, column_gap = gaps,
  cluster_columns = FALSE, show_column_names = FALSE, name = '\nMorphological \nstage\n',
  col = phenotype_cols,
  heatmap_legend_param = list(
    grid_width = unit(legend_width, "cm"), grid_height = unit(legend_height, "cm"),
    labels_gp = gpar(fontsize = cluster_font), title_gp = gpar(fontsize = cluster_font)
  ),
  row_names_gp = gpar(fontsize = cluster_font),
  top_annotation = HeatmapAnnotation(
  # Phenotype group annotation
  Phenotype = anno_block(gp = gpar(fill = group_cols_ordered, col = group_cols_ordered, lwd = 2,
                                   linejoin = 'mitre'),
                         labels = groups_ordered,
                         labels_gp = gpar(col = text_cols_ordered, fontsize = cluster_font),
  ),
  # Strain heatmap annotations
  Strain = anno_block(gp = gpar(fill = 'white', col = group_cols_ordered, lwd = 2,
                                linejoin = 'mitre'),
                      labels = strain_order,
                      labels_gp = gpar(col = 'black', fontsize = cluster_font),
  ),
  annotation_name_gp = gpar(fontsize = cluster_font)
  )
  )
  return(ht_list)
}

# Milestone data
data_impulse <- read.table(paste(path_milestones, 'DEacrossMainStages_AX4_summary_fdr0.001.tsv', sep = ''),
                           header = TRUE, sep = '\t', row.names = 1)
data_deseq <- read.table(paste(path_milestones, 'DE_combined.tsv', sep = ''),
                         header = TRUE, sep = '\t', row.names = 1)
data <- merge(data_deseq, data_impulse, all = TRUE, by = "row.names")
row.names(data) <- data$Row.names
data <- data[, colnames(data) != 'Row.names']

# Expression scale range
expressions <- within(avg_expression_stages, rm('Strain', 'Group', 'main_stage'))
min_expression <- min(expressions)
max_expression <- max(expressions)

# List computed comparisons between stages
comparisons <- c()
for (col in colnames(data)) {
  if (grepl('_FDR_overall', col, fixed = TRUE)) {
    comparison <- gsub('_FDR_overall', '', col, fixed = TRUE)
    comparisons <- append(comparisons, comparison)
  }
}

# Make lists of genes for each comparison and up/down regulation. Order with optical clustering.

milestone_lists <- list()
for (comparison in comparisons) {
  fc_col <- paste(comparison, '_log2FoldChange', sep = '')
  fdr_col <- paste(comparison, '_FDR_overall', sep = '')
  data_defined <- data[!is.na(data[fc_col]) &
                         !is.na(data[fdr_col]) &
                         !is.na(data[comparison]),]
  genes_up <- rownames(data_defined[data_defined[fc_col] >= 2 &
                                      data_defined[fdr_col] <= 0.01 &
                                      data_defined[comparison] == 1,])
  genes_down <- rownames(data_defined[data_defined[fc_col] <= -2 &
                                        data_defined[fdr_col] <= 0.01 &
                                        data_defined[comparison] == 1,])
  genes_up <- optically_order_genes(genes = genes_up, avg_expression = avg_expression_stages)
  genes_down <- optically_order_genes(genes = genes_down, avg_expression = avg_expression_stages)
  milestone_lists[[comparison]][['up']] <- genes_up
  milestone_lists[[comparison]][['down']] <- genes_down
}

# Milestones heatmap annotation
ht_list <- make_anno_mainstage()

# Concatenated heatmaps for each group of milestone gnes
# Some annotations are added only to the first heatmap
first <- TRUE
for (comparison in comparisons) {
  # Comparison genes
  genes_up <- milestone_lists[[comparison]][['up']]
  genes_down <- milestone_lists[[comparison]][['down']]
  split <- c(rep('down', length(genes_down)), rep('up', length(genes_up)))
  genes <- c(genes_down, genes_up)

  # Edit stage names
  data_anno <- comparison
  data_anno <- gsub('no_agg', 'no agg', data_anno)
  data_anno <- gsub('_', ' to ', data_anno)
  data_anno <- gsub('no agg', 'no_agg', data_anno)

  # Annotation for up and downergulated genes
  direction_annotation <- rowAnnotation(Direction = split,
                                        col = list(Direction = c('down' = '#5673e0', 'up' = '#d44e41')),
                                        show_legend = first, annotation_name_side = "top", show_annotation_name = first,
                                        annotation_name_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                                        annotation_legend_param = list(
                                          Direction = list(title = '\nDirection\n',
                                                           grid_width = unit(legend_width, "cm"),
                                                           grid_height = unit(legend_height, "cm"),
                                                           labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                                                           title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily))
                                        )
  )

  # Milestone comparison expression heatmap
  heatmap <- Heatmap(t(avg_expression_stages[, genes]), cluster_columns = FALSE, cluster_rows = FALSE,
                     show_column_names = FALSE, row_split = split,
                     show_row_names = FALSE, col = viridis(256), column_title = NULL,
                     row_title = data_anno,
                     show_heatmap_legend = first, heatmap_legend_param = list(
      title = "\nRelative \nexpression\n",
      at = c(min_expression, round(mean(c(min_expression, max_expression)), 1), max_expression),
      grid_width = unit(legend_width, "cm"), grid_height = unit(legend_height, "cm"),
      labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
      title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
    ),
                     row_title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                     left_annotation = direction_annotation, row_gap = unit(strain_gap, gap_units)
  )
  first <- FALSE
  ht_list <- ht_list %v% heatmap

}
# Plot and save the heatmap
pdf(paste0(path_save, 'expressionHeatmap_milestones.pdf'), width = 35, height = 30)
draw(ht_list, ht_gap = unit(group_gap, gap_units))
graphics.off()
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

path_data <- 'Data/'
path_save <- 'Results/heatmaps/'
path_avg <- 'Results/averaged/'
if (!dir.exists(path_save)) dir.create(path_save, recursive = TRUE)

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
#' @param only_ax4 If TRUE use only avg_expression rows that have "Strain" column value equal to "AX4"
#' @returns Ordered list of genes
optically_order_genes <- function(genes, avg_expression, only_ax4 = TRUE) {
  if (length(genes) > 2) {
    if (only_ax4) {
      expression <- t(avg_expression[avg_expression$Strain == 'AX4', genes])
    }else {
      expression <- t(avg_expression[, genes])
    }
    distances <- dist(expression, method = "Euclidean")
    hc <- hclust(d = distances, method = "ward.D2")
    hc_ordered <- reorder(x = hc, dist = distances)
    genes <- as.dendrogram(hc_ordered) %>% labels
  }
  return(genes)
}

# *** Data common to regulons and disaggregation heatmaps

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

# Data used for phenotype groups and strains annotation
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

# Annotation for regulons heatmap
make_annotation <- function(phenotypes_font = parent.frame()$phenotypes_font,
                            legend_height = parent.frame()$legend_height,
                            legend_width = parent.frame()$legend_width,
                            top_annotation_height = parent.frame()$top_annotation_height,
                            phenotype_annotation_height = parent.frame()$phenotype_annotation_height,
                            cluster_font = parent.frame()$cluster_font) {

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

# Expression range for legend
expressions <- within(avg_expression, rm('Time', 'Strain', 'Group'))
min_expression <- min(expressions)
max_expression <- max(expressions)

# ********************
# *** Regulons heatmap

# *** Rename regulons so that regulon with lower number will be plotted first.

# Needs to be done only once as the regulons file is overwritten
path_regulons <- 'Results/regulons/'
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
  clusters <- unique(regulons$Cluster)
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

for (regulons_file in c("regulonsAX4.tab",
                        'regulonsAll.tab')) {
  # Load regulons
  # Regulon groups tab file: First column lists genes and
  # a column named Cluster specifying cluster/regulon of each gene
  regulons <- read.table(paste(path_regulons, regulons_file, sep = ''), header = TRUE, sep = "\t")
  #Name the first column (should contain genes
  colnames(regulons)[1] <- 'Gene'

  # Rename the regulons in plotting order
  cluster_order <- sort_clusters(regulons = regulons, expression_patterns = expression_patterns, pattern_type = 'Peak')
  # Rename the regulons based on the new ordering and save the renamed regulons
  cluster_map <- c(paste('C', c(1:nrow(cluster_order)), sep = ''))
  names(cluster_map) <- as.vector(cluster_order$Cluster)
  remap_cluster <- function(x) { return(cluster_map[[x]]) }
  regulons['Cluster'] <- unlist(map(as.character(regulons$Cluster), remap_cluster))

  # Save renamed regulons
  write.table(regulons, paste(path_regulons, regulons_file, sep = ''), row.names = FALSE, sep = "\t")
}

# *** Load data common to all regulons heatmaps

# Regulons reference (AX4) groups tab file (used for side annotation): First column lists genes and
# a column named Cluster specifying cluster/regulon of each gene
regulons2 <- read.table(paste(path_regulons,"regulonsAX4.tab",
                              sep = ''), header = TRUE, sep = "\t")
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

# Reference (AX4) regulon colours: 13 distinct colours ordered by rainbow
colours_regulons2 <- c('#800000', '#e6194b', '#f58231', '#9a6324', '#ffe119', '#9dd100', '#3cb44b', '#aaffc3',
                       '#46f0f0', '#2e97ff', '#000075', '#911eb4', '#f032e6','#ed95c5')
# Refernce (AX4) regulons colour map
colours_regulons2_map <- colours_regulons2[1:length(unique(regulons2$Cluster))]
if(length(unique(regulons2$Cluster)) >length(colours_regulons2)){
  print('Number of availiable colours is smaller than number of regulons. Add colours to variable colours_regulons2.')
}
regulons2_clusters <- as.vector(unique(regulons2$Cluster))
names(colours_regulons2_map) <- regulons2_clusters[order(nchar(regulons2_clusters), regulons2_clusters)]

# *** Plot AX4 and all strains based regulons
# Plot all strains based or AX4 based regulons.
regulons_data <- list(all = list(
  file = 'regulonsAll.tab', letters = FALSE),
                      'AX4' = list(file = "regulonsAX4.tab", letters = TRUE))
for (regulons_data_sub in regulons_data) {
  regulons_file <- regulons_data_sub$file
  regulon_to_letters <- regulons_data_sub$letters

  # Regulon groups tab file: First column lists genes and
  # a column named Cluster specifying cluster/regulon of each gene
  regulons <- read.table(paste(path_regulons, regulons_file, sep = ''), header = TRUE, sep = "\t")
  # Name the first column (should contain genes)
  colnames(regulons)[1] <- 'Gene'

  # Get regulons - list unique and sort
  clusters <- unique(regulons$Cluster)
  vals <- as.numeric(gsub("C", "", clusters))
  clusters <- clusters[order(vals)]

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
    # Rename cluster number to a letter (for AX4)
    if (regulon_to_letters) cluster_anno <- LETTERS[as.integer(cluster_anno)]

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
  draw(ht_list)
  graphics.off()
}
# ***************************
# *** Disaggregation heatmaps
path_disaggregation <- 'Results/disaggregation/'

# *** Heatmap for genes DE in tgrB1 and tgrB1C1, but not in AX4, comH, or tagB

# DE filtering thresholds
padj <- 0.01
lfc <- 1.32

ht_list <- make_annotation()
first <- TRUE
for (times in list(c(4, 6), c(6, 8), c(8, 12))) {
  # Find genes DE in tgrB1 and tgrB1C1, but not in AX4, comH, or tagB for each time comparison
  data_AX4 <- read.table(paste0(path_disaggregation, 'AX4_confoundRep_FDRoptim0.01_DE_hr', times[2],
                                '_ref_hr', times[1], '.tsv'),
                         sep = '\t', header = TRUE, row.names = 1)
  data_comH <- read.table(paste0(path_disaggregation, 'comH_confoundRep_FDRoptim0.01_DE_hr', times[2],
                                 '_ref_hr', times[1], '.tsv'),
                          sep = '\t', header = TRUE, row.names = 1)
  data_tagB <- read.table(paste0(path_disaggregation, 'tagB_confoundRep_FDRoptim0.01_DE_hr', times[2],
                                 '_ref_hr', times[1], '.tsv'),
                          sep = '\t', header = TRUE, row.names = 1)
  data_tgrB1 <- read.table(paste0(path_disaggregation, 'tgrB1_confoundRep_FDRoptim0.01_DE_hr', times[2],
                                  '_ref_hr', times[1], '.tsv'),
                           sep = '\t', header = TRUE, row.names = 1)
  data_tgrB1C1 <- read.table(paste0(path_disaggregation, 'tgrB1C1_confoundRep_FDRoptim0.01_DE_hr', times[2],
                                    '_ref_hr', times[1], '.tsv'),
                             sep = '\t', header = TRUE, row.names = 1)

  filter <- Reduce(union, list(
    rownames(data_AX4)[data_AX4$padj <= padj & data_AX4$log2FoldChange >= lfc],
    rownames(data_tagB)[data_tagB$padj <= padj & data_tagB$log2FoldChange >= lfc],
    rownames(data_comH)[data_comH$padj <= padj & data_comH$log2FoldChange >= lfc]))
  genes <- Reduce(intersect, list(
    rownames(data_tgrB1)[data_tgrB1$padj <= padj & data_tgrB1$log2FoldChange >= lfc],
    rownames(data_tgrB1C1)[data_tgrB1C1$padj <= padj & data_tgrB1C1$log2FoldChange >= lfc]
  ))
  genes <- genes[!(genes %in% filter)]
  # Sort genes
  genes <- optically_order_genes(genes = genes, avg_expression = avg_expression)

  heatmap <- Heatmap(t(avg_expression[, genes]), cluster_columns = FALSE, cluster_rows = FALSE,
                     show_column_names = FALSE, show_row_names = FALSE, col = viridis(256), column_title = NULL,
                     row_title = paste0(times[2], 'hr vs ', times[1], 'hr'), show_heatmap_legend = first,
                     heatmap_legend_param = list(
                       title = "\nRelative \nexpression\n",
                       at = c(min_expression, round(mean(c(min_expression, max_expression)), 1), max_expression),
                       grid_width = unit(legend_width, "cm"), grid_height = unit(legend_height, "cm"),
                       labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                       title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)),
                     row_title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),)
  first <- FALSE
  ht_list <- ht_list %v% heatmap
}

#Plots the heatmap
pdf(paste0(path_save,
           'expressionHeatmap_tgrB1tgrB1C1vsAX4tagBcomH_confoundRep_FDRoptim0.01_DEpadj', padj,
           'lfc', lfc, '.pdf'), width = 35, height = 20)
draw(ht_list)
graphics.off()

# *** Heatmap for genes upregulated in tag_dis group compared to whole AX4 profile
ht_list <- make_annotation()

# Find DE genes
data <- read.table(paste0(path_disaggregation,
                          'tagdisVSAX4all_alternativegreater_FDRoptim0.01_DE_tgrB1hr8hr10hr12andtgrB1C1hr8hr10hr12_ref_AX4all.tsv'),
                   sep = '\t', header = TRUE, row.names = 1)
padj <- 0.01
lfc <- 2
genes <- rownames(data)[data$padj <= padj & data$log2FoldChange >= lfc]
# Sort genes
genes <- optically_order_genes(genes = genes, avg_expression = avg_expression)

heatmap <- Heatmap(t(avg_expression[, genes]), cluster_columns = FALSE, cluster_rows = FALSE, show_column_names = FALSE,
                   show_row_names = FALSE, col = viridis(256), column_title = NULL, show_heatmap_legend = TRUE,
                   heatmap_legend_param = list(
                     title = "\nRelative \nexpression\n",
                     at = c(min_expression, round(mean(c(min_expression, max_expression)), 1), max_expression),
                     grid_width = unit(legend_width, "cm"), grid_height = unit(legend_height, "cm"),
                     labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                     title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
                   )
)
ht_list <- ht_list %v% heatmap

# Plot the heatmap
pdf(paste0(path_save,
           'expressionHeatmap_tagdisVSAX4all_alternativegreater_FDRoptim0.01DE_tgrB1hr8hr10hr12andtgrB1C1hr8hr10hr12_ref_AX4all_DEpadj', padj,
           'lfc', lfc, '.pdf'), width = 35, height = 8)
draw(ht_list)
graphics.off()

# *** Heatmap for genes obtained on data from Nichols, et al. (2020)

# *** Find upregulated genes on data from Nichols, et al. (2020)
data_mb <- read.table(paste0(path_disaggregation,
                             'mediaVSbufferAll_alternativegreater_FDRoptim0.01_DE_media0.5hr1hr2hr_ref_bufferAll.tsv'),
                      sep = '\t', header = TRUE, row.names = 1)
padj <- 0.01
lfc <- 2
genes <- rownames(data_mb)[data_mb$padj <= padj & data_mb$log2FoldChange >= lfc]
genes <- optically_order_genes(genes = genes, avg_expression = avg_expression)

# *** Draw heatmap on data used for other heatmaps
ht_list <- make_annotation()
heatmap <- Heatmap(t(avg_expression[, genes]), cluster_columns = FALSE, cluster_rows = FALSE,
                   show_column_names = FALSE, show_row_names = FALSE, col = viridis(256),
                   column_title = NULL, show_heatmap_legend = TRUE,
                   heatmap_legend_param = list(
                     title = "\nRelative \nexpression\n",
                     at = c(min_expression, round(mean(c(min_expression, max_expression)), 1), max_expression),
                     grid_width = unit(legend_width, "cm"), grid_height = unit(legend_height, "cm"),
                     labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                     title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
                   )
)
ht_list <- ht_list %v% heatmap

pdf(paste0(path_save,
           'expressionHeatmap_mediaVSbufferAll_alternativegreater_FDRoptim0.01_DE_media0.5hr1hr2hr_ref_bufferAll_padj',
           padj, 'lfc', lfc, '.pdf'), width = 35, height = 8)
draw(ht_list)
graphics.off()

# *** Draw heatmap on data from Nichols, et al. (2020)

# Load new expression data and reorder genes
avg_expression_mb <- read.table(paste0(path_avg, "genesMediaBuffer_averaged_scaled_percentile99_max0.1.tsv"),
                                header = TRUE, row.names = 1, sep = "\t")
genes <- optically_order_genes(genes = genes, avg_expression = avg_expression_mb, only_ax4 = FALSE)

# Header annotation
# Time annotation
times_mb <- unique(avg_expression_mb['Time'])
col_time_mb <- colorRamp2(c(min(times_mb), max(times_mb)), c("white", "#440154FF"))
ht_list <- Heatmap(t(avg_expression_mb['Time']), height = unit(top_annotation_height, "cm"),
                   # Split by media and buffer
                   column_split = factor(avg_expression_mb$Group), column_title = NULL,
                   cluster_columns = FALSE, cluster_rows = FALSE,
                   show_column_names = FALSE, name = '\nTime\n', col = col_time_mb,
                   heatmap_legend_param = list(
                     at = c(min(times_mb), as.integer(mean(c(min(times_mb), max(times_mb)))), max(times_mb)),
                     grid_width = unit(legend_width, "cm"), grid_height = unit(legend_height, "cm"),
                     labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                     title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
                   ),
                   row_names_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                   # Group annotation
                   top_annotation = HeatmapAnnotation(
                     Group = anno_block(gp =
                                          gpar(
                                            fill = c('#00b2ff', '#ed1c24'), col = c('#00b2ff', '#ed1c24'),
                                            lwd = 2, linejoin = 'mitre'
                                          ),
                                        labels = c('buff', 'media'),
                                        labels_gp = gpar(col = 'black', fontsize = cluster_font, fontfamily = fontfamily),
                                        show_name = TRUE),
                     annotation_name_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
                   )
)

# Expression range for legend
expressions_mb <- within(avg_expression_mb, rm('Time', 'Group'))
min_expression_mb <- min(expressions_mb)
max_expression_mb <- max(expressions_mb)

# Expression heatmap
heatmap <- Heatmap(t(avg_expression_mb[, genes]), cluster_columns = FALSE, cluster_rows = FALSE,
                   show_column_names = FALSE, show_row_names = FALSE, col = viridis(256), column_title = NULL,
                   show_heatmap_legend = TRUE,
                   heatmap_legend_param = list(
                     title = "\nRelative \nexpression\n",
                     at = c(min_expression_mb, round(mean(c(min_expression_mb, max_expression_mb)), 1), max_expression_mb),
                     grid_width = unit(legend_width, "cm"), grid_height = unit(legend_height, "cm"),
                     labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                     title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
                   )
)

# Combine annotation and expression heatmaps
ht_list <- ht_list %v% heatmap
pdf(paste0(path_save,
           'expressionHeatmapMediaBuffer_mediaVSbufferAll_alternativegreater_FDRoptim0.01_DE_media0.5hr1hr2hr_ref_bufferAll_padj',
           padj, 'lfc', lfc, '.pdf'), width = 10, height = 8)
draw(ht_list)
graphics.off()

# ***************************
# *** UNUSED Aberrant neighbourhood heatmaps
#if(FALSE){
# *** Prepare data
path_aberrant <- 'Results/aberrant_neighbourhood/'

# Median similarities to AX4 based neighbours
median_sims <- read.table(paste0(path_aberrant, 'simsMedian_AX4basedNeighByStrain_kN11_samples10resample10.tsv'),
                          header = TRUE, sep = '\t', row.names = 1)
median_sims <- median_sims[, strain_order]

# Ordered strain groups
strain_groups <- c('prec', 'sFB', 'cud', 'tag', 'tag_dis', 'lag_dis', 'agg-')

# Selecte genes with aberrant neighbourhood in strain groups
# FDR from U test
FDR <- 0.001
# Difference between medians of neighbourhood similarities of the two groups
MEDIFF <- 0.3
# Select genes - a DF with groups in columns, genes in row, and values representing selection (1 - yes, 0 - no)
all_abberant <- NULL
for (group in strain_groups) {
  data <- read.table(paste0(path_aberrant, 'comparisonsSims_', group,
                            '_AX4basedNeighByStrain_kN11_samples10resample10.tsv'),
                     header = TRUE, sep = '\t')
  genes <- as.vector(data[data$FDR <= FDR & data['Difference.median'] >= MEDIFF, 'Gene'])
  n_genes <- length(genes)
  df <- data.frame('Gene' = genes, group = rep(1, (length(genes))))
  colnames(df) <- c('Gene', group)
  if (is.null(all_abberant)) {
    all_abberant <- df
  }else {
    all_abberant <- merge(all_abberant, df, by = 'Gene', all = TRUE)
  }

}
all_abberant[is.na(all_abberant)] <- 0
rownames(all_abberant) <- all_abberant$Gene
all_abberant <- all_abberant[, strain_groups]

# List specifiing for each group how to colour unselected (white) and selected genes (strain group colour)
selected_cols <- list()
for (group in strain_groups) {
  selected_cols[[group]] <- c('0' = 'white', '1' = group_cols[[group]])
}
 # Ordered unique group colours
group_cols_ordered_unique<-(unlist(lapply(strain_groups, function(x) group_cols[[x]])))
# *** Expression hetamap of selected genes in each group

for (group in strain_groups) {

  # Select and order genes of a group
  genes <- rownames(all_abberant)[all_abberant[group] == 1]
  #print(paste('Genes with aberrant neighbourhood: N selected genes in', group, ':', length(genes)))
  genes <- optically_order_genes(genes = genes, avg_expression = avg_expression)

  # Annotate rows with information about in which group genes were termed to have aberrant neighbourhood
  row_annotation <- rowAnnotation(df = all_abberant[genes,], col = selected_cols, gp = gpar(col = NA),
                                  show_legend = FALSE, annotation_name_side = "top", show_annotation_name = TRUE,
                                  annotation_name_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
  )

  # Clustom legend for selected genes annotation
  # Displays strain group colours, denoting that a gene was abberantly expressed
  lgd_list <- list(Legend(
    labels = strain_groups, title = "Aberrant\nneighborhood\n",
    legend_gp = gpar(col = group_cols_ordered_unique, fill = group_cols_ordered_unique, fontsize = cluster_font, fontfamily = fontfamily),
    grid_width = unit(legend_width, "cm"), grid_height = unit(legend_height, "cm"),
    labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
    title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
  ))

  # Expression heatmap
  ht_list <- make_annotation()

  heatmap <- Heatmap(t(avg_expression[, genes]), cluster_columns = FALSE, cluster_rows = FALSE,
                     show_column_names = FALSE, show_row_names = FALSE, col = viridis(256), column_title = NULL,
                     show_heatmap_legend = TRUE,
                     heatmap_legend_param = list(
                       title = "\nRelative \nexpression\n",
                       at = c(min_expression, round(mean(c(min_expression, max_expression)), 1), max_expression),
                       grid_width = unit(legend_width, "cm"), grid_height = unit(legend_height, "cm"),
                       labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                       title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
                     ),
                     left_annotation = row_annotation
  )
  ht_list <- ht_list %v% heatmap

  #Plots the combined heatmap and legend
  pdf(paste0(path_save, group, '_padj', FDR, 'MEdiff', MEDIFF,
             '_AX4basedNeighByStrain_kN11_samples10resample10.pdf'),
      width = 40, height = 25)
  draw(ht_list, annotation_legend_list = lgd_list)
  graphics.off()
}

# *** Heatmap of ordered selected genes and median neighbourhood similarities of each strain

# Order genes so that selected gene annotation is optimally ordered
distances <- dist(all_abberant, method = "binary")
hc <- hclust(d = distances, method = "ward.D2")
hc_ordered <- reorder(x = hc, dist = distances)
genes <- as.dendrogram(hc_ordered) %>% labels

# DF with gene selection being marked with colours
# As all_abberant, but with NA instead of 0 and group colour instead of 1
abberant_colours <- all_abberant[genes,]
for (col in colnames(abberant_colours)) {
  new_col <- abberant_colours[col]
  new_col[new_col == '1'] = col
  new_col[new_col == '0'] = NA
  abberant_colours[col] <- new_col
}

# Heatmap indicating which gene had abberant neighbourhood in which strain group
heatmap_selected <- Heatmap(as.matrix(abberant_colours),
                            cluster_columns = FALSE, cluster_rows = FALSE,
                            col = group_cols, show_row_names = FALSE, na_col = 'white',
                            heatmap_legend_param = list(
                              title = "\nAberrant\nneighbourhood\n",
                              grid_width = unit(legend_width, "cm"), grid_height = unit(legend_height, "cm"),
                              labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                              title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
                            ),
                            column_title_gp = (gpar(fontsize = cluster_font, fontfamily = fontfamily)),
                            width = 3, column_names_side = "top"
)


# Heatmap of median similarities to gene neighbours across strains
heatmap_sims <- Heatmap(as.matrix(median_sims[genes,]), col = rev(viridis(256)),
                        column_split = factor(colnames(median_sims), levels = strain_order),
                        column_gap = unit(
                          gsub(gap_units, '',
                               gsub(unit(strain_gap, gap_units), unit(0, gap_units), gaps, fixed = TRUE)
                          ), gap_units),
                        column_title = NULL, cluster_columns = FALSE, cluster_rows = FALSE,
                        show_row_names = FALSE, show_column_names = FALSE,
                        heatmap_legend_param = list(
                          title = "\nMedian\nsimilarity\n",
                          grid_width = unit(legend_width, "cm"), grid_height = unit(legend_height, "cm"),
                          labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
                          title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
                        ),
                        # Annotation with strain groups and strains, as for expression heatmaps
                        # Strains in median similarities DF are ordered as in average expression DF
                        top_annotation = HeatmapAnnotation(
                          Phenotype = anno_block(
                            gp =gpar(fill = group_cols_ordered, col = group_cols_ordered, lwd = 2, linejoin = 'mitre'),
                            labels = groups_ordered,
                            labels_gp = gpar(col = text_cols_ordered, fontsize = cluster_font, fontfamily = fontfamily),
                            show_name = TRUE
                          ),
                          Strain = anno_block(
                            gp =gpar(fill = 'white', col = group_cols_ordered, lwd = 2, linejoin = 'mitre'),
                            labels = strain_order,
                            labels_gp = gpar(col ='black', fontsize = cluster_font, fontfamily = fontfamily ),
                            show_name = TRUE
                          ),
                          annotation_name_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
                        )
)
pdf(paste0(path_save, 'selectedSimilarities_padj', FDR, 'MEdiff', MEDIFF,
           '_AX4basedNeighByStrain_kN11_samples10resample10.pdf'), width = 33, height = 25)
draw(heatmap_selected + heatmap_sims)
graphics.off()

# *** Heatmap of overlap between groups
#' Overlap between two vectors normalised by the length of the smaller vector
#' @param s1 vector 1
#' @param s2 vector2
#' @return len_overlap/len_smaller_vector
proportion_smaller<-function(s1,s2){
    return(length(intersect(s1,s2))/min(length(s1),length(s2)))
}

# Matrix with relative overlap between selected genes in each group
n_groups<-length(strain_groups)
group_overlap<-matrix(,nrow=n_groups,ncol=n_groups)
for(i in 1:(n_groups-1)){
    for(j in (1+i):n_groups){
      genes1 <- rownames(all_abberant)[all_abberant[strain_groups[i]] == 1]
      genes2<-rownames(all_abberant)[all_abberant[strain_groups[j]] == 1]
      group_overlap[i,j]=proportion_smaller(genes1,genes2)
    }
}
colnames(group_overlap)<-strain_groups
rownames(group_overlap)<-strain_groups

overlap_heatmap<-Heatmap(group_overlap,cluster_columns = FALSE,cluster_rows = FALSE,
        col = viridis(256),
        heatmap_legend_param = list(
          title = "Ratio\nsmaller\noverlap\n",
          grid_width = unit(legend_width, "cm"), grid_height = unit(legend_height, "cm"),
          labels_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily),
          title_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily)
        ),
        column_names_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily, col=group_cols_ordered_unique ),
        row_names_gp = gpar(fontsize = cluster_font, fontfamily = fontfamily,col=group_cols_ordered_unique )
)
pdf(paste0(path_save, 'selectedOverlap_padj', FDR, 'MEdiff', MEDIFF,
           '_AX4basedNeighByStrain_kN11_samples10resample10.pdf'), width = 8, height = 7)
draw(overlap_heatmap)
graphics.off()
#}
# ***************************
# ****** Milestones heatmap

# *** Load data and helper functions
path_milestones <- 'Results/milestones/'

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
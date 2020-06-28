print('Preparing disaggregation data.')

source('helper_DE.R')

# *****************
# *** Load data

path_data <- 'Data/'
path_save <- 'Results/disaggregation/'
if (!dir.exists(path_save)) dir.create(path_save, recursive = TRUE)

genes <- read.table(paste(path_data, "mergedGenes_counts.tsv", sep = ''), header = TRUE, row.names = 1, sep = "\t")
conditions <- read.table(paste(path_data, "conditions_mergedGenes.tsv", sep = ''), header = TRUE,
                         row.names = 'Measurment', sep = "\t")
# Gene data frame column names (sample names) were modified during import, so chnage the metadata accordingly
rownames(conditions) <- make.names(rownames(conditions))

# **********************
# *** DE between timepoints of a strain
for (strain in c('tgrB1', 'tgrB1C1', 'AX4', 'comH', 'tagB')) {
  for (timepoints in list(c(4, 6), c(6, 8), c(8, 12))) {

    # Get samples of both timepoints and prepare data for DESeq2
    print(paste('Analysisng', paste0(strain, ' ', timepoints, 'hr', collapse = ' vs ')))
    name1 <- paste0('hr', timepoints[1])
    name2 <- paste0('hr', timepoints[2])
    samples1 <- conditions[conditions$Strain == strain & conditions$Time == timepoints[1], c('Replicate', 'Time')]
    samples1['Comparison'] <- name1
    samples2 <- conditions[conditions$Strain == strain & conditions$Time == timepoints[2], c('Replicate', 'Time')]
    samples2['Comparison'] <- name2
    # Keep only replicates present in both groups/timepoints
    replicates_both <- intersect(unique(samples1$Replicate), unique(samples2$Replicate))
    conditions_sub <- rbind(samples1, samples2)
    conditions_sub <- conditions_sub[conditions_sub$Replicate %in% replicates_both,]
    print(paste('N samples 1', sum(conditions_sub$Comparison == name1),
                'N samples 2', sum(conditions_sub$Comparison == name2)))
    genes_sub <- genes[, rownames(conditions_sub)]

    # Make DESeq2 comparison
    # DESeq2 design
    design_formula <- ~Replicate + Comparison
    # Adjust for the used padj threshold
    fdr_optim <- 0.01
    path_save_comparison <- paste0(path_save, strain, '_confoundRep_FDRoptim', fdr_optim, '_')
    runDeSeq2(conditions = conditions_sub, genes = genes_sub, case = name2, control = name1,
              design = design_formula, main_lvl = 'Comparison',
              path = path_save_comparison, alpha = fdr_optim
    )
  }
}

# **********************
# *** DE between specific tag_dis timepoints and all AX4 timepoints

# Prepare data for DESeq2
timepoints <- c(8, 10, 12)
times_str <- paste0('hr', timepoints, collapse = '')
name1 <- paste0('tgrB1', times_str, 'andtgrB1C', times_str)
name2 <- "AX4all"
print(paste('Comparing', name2, 'vs', name1))
samples1 <- conditions[(conditions$Strain == 'tgrB1' & conditions$Time %in% timepoints) |
                         (conditions$Strain == 'tgrB1C1' & conditions$Time %in% timepoints), c('Strain', 'Time')]
samples2 <- conditions[conditions$Strain == 'AX4', c('Strain', 'Time')]
samples1['Comparison'] <- name1
samples2['Comparison'] <- name2
print(paste('N samples 1', nrow(samples1), 'N sam,ples 2', nrow(samples2)))
conditions_sub <- rbind(samples1, samples2)
genes_sub <- genes[, rownames(conditions_sub)]

# Run DESeq2
design_formula <- ~Comparison
# Alternative hypothesis and FDR optimisation
alternative <- 'greater'
fdr_optim <- 0.01
path_save_comparison <- paste0(path_save, 'tagdisVSAX4all_alternative', alternative, '_FDRoptim', fdr_optim, '_')
runDeSeq2(conditions = conditions_sub, genes = genes_sub, case = name1, control = name2,
          design = design_formula, main_lvl = 'Comparison', path = path_save_comparison,
          altHypothesis = alternative, alpha = fdr_optim
)

# **********************
# *** DE genes in data from Nichols, et al. (2020)
print('Performing analysis for data from Nichols, et al. (2020)')
# *** Load data from Nichols, et al. (2020)
genes_mb <- read.table(paste(path_data, "mediaBuffer_counts.tsv", sep = ''), header = TRUE, row.names = 1, sep = "\t")
conditions_mb <- read.table(paste(path_data, "conditions_mediaBuffer.tsv", sep = ''), header = TRUE,
                            row.names = 'Measurment', sep = "\t")
rownames(conditions_mb) <- make.names(rownames(conditions_mb))

# *** DE analysis - compare specific media timepoints towards all buffer samples
# Media timepoints
times <- c(0.5, 1, 2)
# Prepare DESeq2 input data
name1 <- paste0('media', paste(times, collapse = 'hr'), 'hr')
name2 <- "bufferAll"
samples1 <- conditions_mb[(conditions_mb$Group == 'media' & conditions_mb$Time %in% times), c('Group', 'Time')]
samples2 <- conditions_mb[conditions_mb$Group == 'buff', c('Group', 'Time')]
samples1['Comparison'] <- name1
samples2['Comparison'] <- name2
print(paste('N samples 1', nrow(samples1), 'N sam,ples 2', nrow(samples2)))
conditions_sub <- rbind(samples1, samples2)
genes_sub <- genes_mb[, rownames(conditions_sub)]
design_formula <- ~Comparison
# Alternative hypothesis and FDR optimisation
alternative <- 'greater'
fdr_optim <- 0.01

path_save_comparison <- paste0(path_save, 'mediaVSbufferAll_',
                               'alternative', alternative, '_FDRoptim', fdr_optim, '_')
runDeSeq2(conditions = conditions_sub, genes = genes_sub, case = name1, control = name2,
          design = design_formula, main_lvl = 'Comparison', path = path_save_comparison,
          altHypothesis = alternative, alpha = fdr_optim
)

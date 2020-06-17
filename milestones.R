print('Preparing milestones data.')

library("DESeq2")
library('dplyr')
library("BiocParallel")
library("ImpulseDE2")
library('plyr')
library(parallel)
# ****************************
# *** Helper functions for DE analysis

  #' Builds a DDS for DESeq2
  #' @param conditions (object class data.frame) Metainformation about expression samples.
  #' Data frame of size M*D, where M are samples/measurements and D are description columns.
  #' @param genes (object class data.frame) Expression data. Data frame of size G*M, where G are genes and M are
  #' samples/measurements.
  #' @param case Value of main_lvl column indicating that sample should be used as case.
  #' @param ref Value of main_lvl column indicating that sample should be used as reference.
  #' @param design (formula) DeSeq2 design formula.
  #' @param main_lvl Column to be used for distinguishing between case and reference samples.
  #' @param filter (int) Remove genes with less than that many counts across all samples. If NULL no
  #' filtering is performed.
  #' @returns (list) Elements: dds - DESeq2 dds DESeqDataSet, coldata - coldata as prepared for DESeq2.
buildDDS <- function(conditions, genes, case, ref, design, main_lvl, filter = 1) {

  genes <- genes[, unlist(conditions[main_lvl], use.names = FALSE) %in% c(case, ref)]
  conditions <- conditions[unlist(conditions[main_lvl], use.names = FALSE) %in% c(case, ref),]

  coldata <- data.frame(lapply(conditions, as.factor))
  coldata <- data.frame(lapply(coldata, droplevels))
  rownames(coldata) <- rownames(conditions)

  #Check that sample descriptions match betweeen expression  data and description
  if (all(rownames(coldata) == colnames(genes))) {
    #Make DeSeq2 data object
    dds <- DESeqDataSetFromMatrix(countData = genes,
                                  colData = coldata,
                                  design = design)
    dds[[main_lvl]] <- relevel(dds[[main_lvl]], ref = ref)
    if (!is.null(filter)) {
      keep <- rowSums(counts(dds)) >= filter
      dds <- dds[keep,]
    }
    return(list(dds = dds, coldata = coldata))
  }
  else {
    stop("Sample names in conditions and genes do not match")
  }
}

#' Test DE with DESeq2, filter, and save or return results.
#' Results dataframe is filtered to remove NA entries and ordered by ascending padj.
#' @param dds (object class DESeqDataSet) DESeq2 dds for analysis.
#' @param sample Case level used in the DE comparison.
#' @param ref Reference level used in the DE comparison.
#' @param path (str) Directory path for saving the results, ending with '/'. The file is saved within the directory
#' with name DE_sample_ref_ref.tsv, where sample and second ref are parameters.
#' If NULL results are not saved and are instead returned.
#' @param removeNA Remove any rows in results that contain NA.
#' @returns Processed results dataframe; if path is NULL.
testDE <- function(dds, sample, ref, main_lvl, path = NULL, removeNA = TRUE) {
  res <- results(dds, contrast = c(main_lvl, sample, ref), parallel = TRUE)
  print(summary(res))
  if (removeNA) {
    resNN <- na.omit(res)
  }else {
    resNN <- res
  }
  resOrder <- resNN[order(resNN$padj),]
  if (is.null(path)) {
    return(resOrder)
  }else {
    write.table(x = resOrder, file = paste(path, 'DE_', sample, '_ref_', ref, '.tsv', sep = ''),
                sep = '\t', col.names = NA)
  }
}

#' Builds DESeq2 object and performs DE analysis based on specified contrast.
#' @param conditions (object class data.frame) Metainformation about expression samples.
#' Data frame of size M*D, where M are samples/measurements and D are description columns.
#' @param genes (object class data.frame) Expression data. Data frame of size G*M, where G are genes and M are
#' samples/measurements.
#' @param case Value of main_lvl column indicating that sample should be used as case.
#' @param control Value of main_lvl column indicating that sample should be used as reference.
#' @param design (formula) DeSeq2 design formula.
#' @param main_lvl Column to be used for distinguishing between case and reference samples.
#' @param path (str) Path for saving the results of DE analysis, as specified in testDE.
#' Name is modified based on testDE function.
#' @param save_dds_path (str) Path for saving DESeqDataSet object after calling DESeq on it. Should be a directory path
#' @param removeNA Remove any rows in results that contain NA.
#' ending in '|'. File is named case_ref_control.rds, where case and control are parameters.
runDeSeq2 <- function(conditions, genes, case, control, design, main_lvl, path = NULL, save_dds_path = NULL,
                      removeNA = TRUE) {
  dds <- buildDDS(conditions = conditions, genes = genes, case = case, ref = control, design = design, main_lvl = main_lvl, filter = 1)$dds
  dds <- DESeq(dds, parallel = TRUE)
  if (!is.null(save_dds_path)) {
    saveRDS(object = dds, file = paste(path, case, '_ref_', control, '.rds', sep = ''))
  }
  if (is.null(path)) return(testDE(dds = dds, sample = case, ref = control, path = path, main_lvl = main_lvl,
                                   removeNA = removeNA))
  else testDE(dds = dds, sample = case, ref = control, path = path, main_lvl = main_lvl, removeNA = removeNA)
}

# *****************
# *** Load data

# Threads
THREADS <- detectCores() - 1
if (THREADS > 30) THREADS <- 30
register(MulticoreParam(THREADS))

path_data <- 'Data/'
path_save <- 'Results/milestones/'
if (!dir.exists(path_save)) dir.create(path_save, recursive = TRUE)

genes <- read.table(paste(path_data, "mergedGenes_counts.tsv", sep = ''), header = TRUE, row.names = 1, sep = "\t")
conditions <- read.table(paste(path_data, "conditions_mergedGenes.tsv", sep = ''), header = TRUE,
                         row.names = 'Measurment', sep = "\t")
# Gene data frame column names (sample names) were modified during import, so chnage the metadata accordingly
rownames(conditions) <- make.names(rownames(conditions))

# Ordered developmental stages
STAGES <- c('no_agg', 'stream', 'lag', 'tag', 'tip', 'slug', 'mhat', 'cul', 'FB')
STAGES_X <- data.frame(Phenotype = STAGES, X = c(1:length(STAGES)))
# ********************
# ******** DE analysis between neighbouring stages (DESeq2)

# Use only AX4 data
conditions_test <- conditions[conditions$Strain == 'AX4',]

# For each pair of neighbouring stages test DE
for (idx in 1:(length(STAGES) - 1)) {
  stage1 <- STAGES[idx]
  stage2 <- STAGES[idx + 1]
  test <- conditions_test[conditions_test['main_stage'] == stage2,]
  control <- conditions_test[conditions_test['main_stage'] == stage1,]
  print(paste('Samples in', stage1, nrow(control), '(control) and in', stage2, nrow(test), '(test)'))

  # Make sure that at least 2 replicates are in control and test group
  if (length(unique(test$Replicate)) > 1 & length(unique(control$Replicate)) > 1) {

    print(paste('Comparing', stage2, 'to', stage1))
    conditions_sub <- rbind(test, control)
    genes_sub <- genes[, rownames(conditions_sub)]
    design_formula <- ~main_stage

    res <- runDeSeq2(conditions = conditions_sub, genes = genes_sub, case = stage2, control = stage1,
                     design = design_formula, main_lvl = 'main_stage', path = path_save, save_dds_path = path_save,
                     removeNA = FALSE)
  }
}

#************************
#******** Identification of gene expression transition times (ImpulseDE2)

# *** Run ImpulseDE2

# Prepare data for ImpulseDE2
Y <- conditions[conditions[, 'main_stage'] != '' & conditions$Strain == 'AX4', 'main_stage', drop = F]
Y[, 'Sample'] <- rownames(Y)
Y[, 'Time'] <- as.numeric(mapvalues(as.character(unlist(Y[, 'main_stage'])), from = STAGES_X$Phenotype, to = STAGES_X$X))
Y[, 'Condition'] <- 'case'
Y <- Y[order(Y$Time),]
X <- genes[, Y$Sample]
print(paste('N samples:', nrow(Y)))

# Run ImpulseDE2
objectImpulseDE2 <- runImpulseDE2(matCountData = as.matrix(X), dfAnnotation = Y, boolCaseCtrl = FALSE,
                                  vecConfounders = NULL, boolIdentifyTransients = TRUE, scaNProc = THREADS)

saveRDS(object = objectImpulseDE2, file = paste(path_save, 'DEacrossMainStages_AX4.rds', sep = ''))
write.table(STAGES_X, file = paste(path_save, 'DEacrossMainStages_AX4_stageOrder.tsv', sep = ''), sep = '\t', row.names = FALSE)

# *** Parse the results of ImpulseDE2 to find expression transitions

# Load ImpulseDE2 and related data
objectImpulseDE2 <- readRDS(paste(path_save, 'DEacrossMainStages_AX4.rds', sep = ''))
stages_x <- read.table(paste(path_save, 'DEacrossMainStages_AX4_stageOrder.tsv', sep = ''), header = TRUE)
stages_x <- stages_x[order(stages_x$X),]
stages_vec <- as.vector(stages_x[, 'Phenotype'])
min_x <- min(stages_x$X)
max_x <- max(stages_x$X)
times_all <- sort(stages_x$X)

# Set ImpulseDE2 DE analysis FDR threshold, default is 0.001
fdr_thershold <- 0.001
objectImpulseDE2 <- updateDEAnalysis(objectImpulseDE2, scaQThresTransients = fdr_thershold)

# Retain only genes defined as DE across stages by ImpulseDE2
model <- objectImpulseDE2@lsModelFits$case
result <- objectImpulseDE2$dfImpulseDE2Results
result <- result[!is.na(result$padj),]
result <- result[result$isMonotonous | result$isTransient,]
annotation <- objectImpulseDE2@dfAnnotationProc
print(paste('N genes DE across stages:', nrow(result)))

# Extrat tranistion times (neighbouring stages)
# Prepare data frame with the following columns: Type - transient or monotnous; Assignmnet - how was the model parsed
# to extratc the transition information (described in more detail below); transition columns of format stage1_stage2
# where 1 indicates transition and 0 no transition; padj values from ImpulseDE2 result: impulse_padj,
# sigmoid_padj (both compared to constant model), impulseTOsigmoid_padj.
parsed <- data.frame()
for (gene in rownames(result)) {
  parsed[gene, c('impulse_padj', 'impulseTOsigmoid_padj', 'sigmoid_padj')] <- result[gene,
                                                                                     c('padj', 'impulseTOsigmoid_padj', 'sigmoidTOconst_padj')]
  type <- ''
  # Was gene determined to be transient or monotonous
  if (result[gene, 'isTransient']) type <- 'transient'
  if (result[gene, 'isMonotonous']) type <- 'monotonous'
  parsed[gene, 'Type'] <- type

  # Parse genes monotonously DE across stages
  if (type == 'monotonous') {

    # If monotonous it can be sigmoid or impulse monotonous. Extract transition time from the model
    # that was better fit compared to the constant model (alternatively, log likelyhood or
    # impulseTOsigmoid_padj could be used).
    # If impulse model was selected the t (transition time) that is closer to the t from monotonous model is used.

    if (result[gene, 'padj'] > result[gene, 'sigmoidTOconst_padj']) {
      model_data <- as.list(model[gene][[1]]$lsSigmoidFit$vecSigmoidParam)
      t <- model_data$t
    }else {
      model_data <- as.list(model[gene][[1]]$lsSigmoidFit$vecSigmoidParam)
      model_impulse <- as.list(model[gene][[1]]$lsImpulseFit$vecImpulseParam)
      diff_t1 <- abs(model_data$t - model_impulse$t1)
      diff_t2 <- abs(model_data$t - model_impulse$t2)
      if (diff_t1 < diff_t2) { t <- model_impulse$t1
      } else { t <- model_impulse$t2 }
    }

    # If t is in between two stages then Assignmnet=impulse.
    # If t is larger/smaller than the last/first stage, respectively, assume that only the last/first stage
    # is on a different expression level (transition is after the first or before the last stage) and
    # Assignment=impulse_manual_border as this is not directly based on the model.
    assignment <- NA
    if (min_x < t & t < max_x) {
      assignment <- 'impulse'
      # Assumes that stages are integers
    } else if (min_x >= t) {
      t <- min_x + 0.1
      assignment <- 'impulse_manual_border'
    } else if (max_x <= t) {
      t <- max_x - 0.1
      assignment <- 'impulse_manual_border'
    }
    parsed[gene, 'Assignment'] <- assignment

    # Between which stages does the transition occur
    phenotypes1 <- as.vector(stages_x[stages_x$X < t, 'Phenotype'])
    phenotypes2 <- as.vector(stages_x[stages_x$X > t, 'Phenotype'])
    transition <- paste(tail(phenotypes1, n = 1), head(phenotypes2, n = 1), sep = '_')
    for (idx in 1:(length(stages_vec) - 1)) {
      neighbouring <- paste(stages_vec[idx], stages_vec[idx + 1], sep = '_')
      if (neighbouring == transition) {
        parsed[gene, neighbouring] <- 1
      }else {
        parsed[gene, neighbouring] <- 0
      }
    }

    # Parse genes transiently DE across stages
  } else if (type == 'transient') {
    model_data <- as.list(model[gene][[1]]$lsImpulseFit$vecImpulseParam)
    t1 <- model_data$t1
    t2 <- model_data$t2
    # t1,t2 can be between different stages - these are then regarded as transition times and Assignment=impulse.
    # Otherwise the model must be further parsed in the order described below.
    # If t1 and t2 are between the same two stages assume that there is a peak or a vallley adjecently to the t1/t2.
    # t1 and t2 are set to be before the first and after the second neighbouring stage, respectively,
    # and Assigmnet=impulse_manual_peak. t1,t2 are also corrected as in impulse_manual_border, if needed.
    # If t1>t2 (except if it was corrected before e.g. in  impulse_manual_peak) t1 and t2 are swapped and
    # Assignment=impulse_manual_swap. t1,t2 are also corrected as in impulse_manual_border, if needed.
    # If t1 or t2 is smaller or larger than the first or last stage, respectively, assume that the transition
    # occurs after the first or before the last stage and Assignment=impulse_manual_border. This can result in t1,t2
    # again being between the same two stages. In this case the model is treated as monotnous and
    # Type=transient_monotonous.
    assignment <- NA
    # t1,t2 are in the right orientation and in the right range for transitions to be determined directly.
    if (min_x < t1 &
      t1 < max_x &
      min_x < t2 &
      t2 < max_x &
      floor(t1) != floor(t2)) {
      assignment <- 'impulse'
    }
    # Check if t1,t2 are between the same stages so that reasigning border does not cause this situation
    if (floor(t1) == floor(t2)) {
      assignment <- 'impulse_manual_peak'

      # Assign t1,t2 so that the two points closest to t1,t2 represent a peak/valley. If this forces
      # t1,t2  below/above min/max stages, respectively, correct as for impulse_manual_border.
      # Assumes that stages/times are integers
      t1 <- max(times_all[times_all < t1]) - 0.1
      t2 <- min(times_all[times_all > t2]) + 0.1
      if (min_x >= t1 & t1 < max_x) {
        t1 <- min_x + 0.1
      }
      if (min_x < t2 & t2 >= max_x) {
        t2 <- max_x - 0.1
      }
    }
    # Check if t2.t1 ans swap t1,t2.
    # If this forces t1,t2 below/above min/max stages, respectively, correct as for impulse_manual_border.
    if (t2 < t1) {
      t_temp <- t1
      t1 <- t2
      t2 <- t_temp
      if (min_x >= t1 & t1 < max_x) {
        t1 <- min_x + 0.1
      }
      if (min_x < t2 & t2 >= max_x) {
        t2 <- max_x - 0.1
      }
      assignment <- 'impulse_manual_swapped'
    }
    # If t1/t2 are below/above min/max stages, respectively, correct as in monotonous model.
    # This may cause t1 and t2 to be equal
    if (min_x >= t1 & t1 < max_x) {
      if (is.na(assignment)) assignment <- 'impulse_manual_border'
      t1 <- min_x + 0.1
    }
    if (min_x < t2 & t2 >= max_x) {
      if (is.na(assignment)) assignment <- 'impulse_manual_border'
      t2 <- max_x - 0.1
    }
    # Check again if t1 and t2 are between the same two stages - the case where this was so originally so was already
    # dealt with. Treat this as monotonous.
    if (floor(t1) == floor(t2)) {
      parsed[gene, 'Type'] <- 'transient_monotonous'
      t2 <- NA
    }

    parsed[gene, 'Assignment'] <- assignment

    # Between which stages does the transition occur
    if (!is.na(assignment)) {
      # Model was parsed as having a single transition
      if (is.na(t2)) {
        phenotypes1 <- as.vector(stages_x[stages_x$X < t1, 'Phenotype'])
        phenotypes2 <- as.vector(stages_x[stages_x$X > t1, 'Phenotype'])
        transition <- paste(tail(phenotypes1, n = 1), head(phenotypes2, n = 1), sep = '_')
        for (idx in 1:(length(stages_vec) - 1)) {
          neighbouring <- paste(stages_vec[idx], stages_vec[idx + 1], sep = '_')
          if (neighbouring == transition) {
            parsed[gene, neighbouring] <- 1
          }else {
            parsed[gene, neighbouring] <- 0
          }
        }
        # Model was parsed as having two transitions
      } else {
        phenotypes1 <- as.vector(stages_x[stages_x$X < t1 | stages_x$X > t2, 'Phenotype'])
        phenotypes2 <- as.vector(stages_x[stages_x$X > t1 & stages_x$X < t2, 'Phenotype'])
        transition_stage1 <- tail(as.vector(stages_x[stages_x$X < t1, 'Phenotype']), n = 1)
        transition_stage2 <- head(as.vector(stages_x[stages_x$X > t1, 'Phenotype']), n = 1)
        transition_stage3 <- tail(as.vector(stages_x[stages_x$X < t2, 'Phenotype']), n = 1)
        transition_stage4 <- head(as.vector(stages_x[stages_x$X > t2, 'Phenotype']), n = 1)
        transition1 <- paste(transition_stage1, transition_stage2, sep = '_')
        transition2 <- paste(transition_stage3, transition_stage4, sep = '_')
        for (idx in 1:(length(stages_vec) - 1)) {
          neighbouring <- paste(stages_vec[idx], stages_vec[idx + 1], sep = '_')
          if (neighbouring == transition1 | neighbouring == transition2) {
            parsed[gene, neighbouring] <- 1
          }else {
            parsed[gene, neighbouring] <- 0
          }
        }
      }
    }
  }
}

# Report on transition parsing
parsed_nonna <- parsed[!is.na(parsed$Assignment),]
print(paste('Assignment: impulse:', sum(parsed_nonna$Assignment == 'impulse'),
            'manual border:', sum(parsed_nonna$Assignment == 'impulse_manual_border'),
            'manual peak:', sum(parsed_nonna$Assignment == 'impulse_manual_peak'),
            'manual swapped:', sum(parsed_nonna$Assignment == 'impulse_manual_swapped'),
            'not-parsed:', sum(is.na(parsed$Assignment))))

# Save the transition information
write.table(parsed, paste(path_save, 'DEacrossMainStages_AX4_summary_fdr', fdr_thershold, '.tsv', sep = ''),
            sep = '\t', col.names = NA)



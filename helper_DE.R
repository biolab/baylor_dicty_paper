library("DESeq2")
library("BiocParallel")
library('dplyr')

# Threads
THREADS <- detectCores() - 1
if (THREADS > 30) THREADS <- 30
register(MulticoreParam(THREADS))

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
#' @param ... Passed to DESeq's results function
testDE <- function(dds, sample, ref, main_lvl, path, removeNA = TRUE, ...) {
  res <- results(dds, contrast = c(main_lvl, sample, ref), parallel = TRUE, ...)
  print(summary(res))
  if (removeNA) {
    resNN <- na.omit(res)
  }else {
    resNN <- res
  }
  resOrder <- resNN[order(resNN$padj),]
  write.table(x = resOrder, file = paste0(path, 'DE_', sample, '_ref_', ref, '.tsv'),
                sep = '\t', col.names = NA)
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
#' @param path (str) Path for saving the results of DE analysis (as specified in testDE) and DESeq object.
#' Name is modified based on case and control name.
#' @param removeNA Remove any rows in results that contain NA.
#' ending in '|'. File is named case_ref_control.rds, where case and control are parameters.
#' @param ... Passed to DESeq's results function
runDeSeq2 <- function(conditions, genes, case, control, design, main_lvl, path , removeNA = TRUE, ...) {
  dds <- buildDDS(conditions = conditions, genes = genes, case = case, ref = control, design = design, main_lvl = main_lvl, filter = 1)$dds
  dds <- DESeq(dds, parallel = TRUE)
  saveRDS(object = dds, file = paste0(path, case, '_ref_', control, '.rds'))
  testDE(dds = dds, sample = case, ref = control, path = path, main_lvl = main_lvl, removeNA = removeNA, ...)
}



library("DESeq2")
library('dplyr')
library("BiocParallel")

# ****************************
# *** Helper functions

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
buildDDS<-function(conditions,genes,case,ref,design,main_lvl,filter=1){

  genes <- genes[, unlist(conditions[main_lvl] , use.names=FALSE) %in% c(case, ref)]
  conditions <- conditions[unlist(conditions[main_lvl]  , use.names=FALSE) %in% c(case, ref),]

  coldata<-data.frame(lapply(conditions, as.factor))
  coldata<-data.frame(lapply(coldata,droplevels))
  rownames(coldata)<-rownames(conditions)

  #Check that sample descriptions match betweeen expression  data and description
  if (all(rownames(coldata)==colnames(genes))){
    #Make DeSeq2 data object
    dds <- DESeqDataSetFromMatrix(countData = genes,
                                  colData = coldata,
                                  design = design)
    dds[[main_lvl]]<- relevel(dds[[main_lvl]], ref = ref)
    if (! is.null(filter)){
      keep <- rowSums(counts(dds)) >= filter
      dds <- dds[keep,]
    }
    return(list(dds=dds,coldata=coldata))
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
#' @returns Processed results dataframe; if path is NULL.
testDE<-function(dds,sample,ref,main_lvl,path=NULL){
  res <- results(dds,contrast=c(main_lvl, sample, ref),parallel = TRUE)
  print(summary(res))
  resNN <- na.omit(res)
  resOrder<-resNN[order(resNN$padj),]
  if (is.null(path)){
    return(resOrder)
  }else {
    write.table(x = resOrder,file =paste( path,'DE_',sample,'_ref_',ref,'.tsv',sep=''),
                sep='\t', col.names=NA)
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
#' ending in '|'. File is named case_ref_control.rds, where case and control are parameters.
runDeSeq2<-function(conditions,genes,case,control,design,main_lvl,path=NULL, save_dds_path=NULL){
  dds<-buildDDS(conditions=conditions,genes=genes,case=case,ref=control,design=design,main_lvl=main_lvl,filter=1)$dds
  dds <- DESeq(dds,parallel = TRUE)
  if(!is.null(save_dds_path)){
      saveRDS(object=dds,file=paste(path,case,'_ref_',control,'.rds',sep=''))
  }
  if (is.null(path)) return(testDE(dds=dds,sample=case,ref=control,path=path,main_lvl=main_lvl))
  else testDE(dds=dds,sample=case,ref=control,path=path,main_lvl=main_lvl)
}

# *****************
# *** Load data
# Threads
register(MulticoreParam(10))

#TODO make data paths the same
path_data1 = '/home/karin/Documents/timeTrajectories/data/countsRaw/combined/'
path_data2 = '/home/karin/Documents/timeTrajectories/data/RPKUM/combined/'
path_save = '/home/karin/Documents/timeTrajectories/data/deTime/neighbouring/AX4/'

genes <- read.table(paste(path_data1, "mergedGenes_counts.tsv", sep = ''), header = TRUE, row.names = 1, sep = "\t")
conditions <- read.table(paste(path_data2, "conditions_mergedGenes.tsv", sep = ''), header = TRUE, row.names = 'Measurment', sep = "\t")
# Gene data frame column names (sample names) were modified during import, so chnage the metadata accordingly
rownames(conditions) <- make.names(rownames(conditions))

# ********************
# ******** DE analysis between neighbouring stages
#Stages are already ordered
STAGES <- c('no_agg', 'stream', 'lag', 'tag', 'tip', 'slug', 'mhat', 'cul', 'FB', 'yem')
PHENOTYPES_ORDERED = c('no_agg', 'stream', 'lag', 'tag', 'tip', 'slug', 'mhat', 'cul', 'FB', 'yem')
PHENOTYPES_X = data.frame(Phenotype = PHENOTYPES_ORDERED, X = c(1:length(PHENOTYPES_ORDERED)))

#Use only AX4 data
conditions_test = conditions[conditions$Strain == 'AX4',]
if (server) {
  path_save = '/home/khrovatin/timeTrajectoriesNet/data/stages/neighbouring/AX4/'
}else {
  path_save = '/home/karin/Documents/timeTrajectories/data/deTime/neighbouring/AX4/'
}

stages = STAGES[STAGES != "yem"]
for (idx in 1:(length(stages) - 1)) {
  stage1 = stages[idx]
  stage2 = stages[idx + 1]
  print(paste('Comparing', stage2, 'to', stage1))
  test <- conditions_test[conditions_test['main_stage'] == stage2,]
  control <- conditions_test[conditions_test['main_stage'] == stage1,]

  #At least 2 replicates are in control and test

  if (length(unique(test$Replicate)) > 1 & length(unique(control$Replicate)) > 1) {

    conditions_sub <- rbind(test, control)
    genes_sub <- genes[, rownames(conditions_sub)]
    design_formula = ~main_stage
    #print(conditions_sub[,c('Replicate','main_stage')])

    res <- runDeSeq2(conditions_sub, genes_sub, case = stage2, control = stage1, design = design_formula, main_lvl = 'main_stage',
                     padj = NULL, logFC = NULL,
                     path = path_save, save_dds_path = path_save)
  }
}
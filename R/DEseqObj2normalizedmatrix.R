## Normalize a RAW COUNTS DEseq2 object (DESeqDataSet) with vst, return the normalized counts matrix.
## Extra parameters (...) are passe to DESeq2::vst()
## Requires DESeq2 and clusterProfiler packages, and functions from customscripts/R/diffexp2gsea.R
source('/home/job/gits/customscripts/R/diffexp2gsea.R')
DE2obj.2.norm.mat <- function(DE2obj = NULL, feature.type = 'ENSEMBL', feature.out = 'SYMBOL', species = 'Homo sapiens', out.dir = getwd(), ...) {
  ## Normalizing (vst)
  DE2obj.norm <- DESeq2::vst(object = DE2obj, ...)
  ## Extracting the normalized count matrix
  norm.mat <- SummarizedExperiment::assay(DE2obj.norm)
  ## Converting feature names if requested
  if(feature.out != feature.type) {
    ## Creating the converting object
    gconv <- as.list(clusterProfiler::bitr(rownames(norm.mat), fromType = feature.type, toType = feature.out, OrgDb=paste0(msigdbr2org(species), '.db')))
    ## Converting features
    names(gconv[[feature.out]]) <- gconv[[feature.type]]
    names(gconv[[feature.type]]) <- gconv[[feature.out]]
    rownames(norm.mat) <- unname(gconv[[feature.out]][as.character(rownames(norm.mat))])
  }
  ## Filtering
  na.tf <- is.na(rownames(norm.mat))
  if(any(na.tf)) norm.mat <- norm.mat[!na.tf,]
  ## Removing duplicates
  dup.tf <- duplicated(rownames(norm.mat))
  if(any(dup.tf)) norm.mat <- norm.mat[!dup.tf,]
  ## Sorting features
  norm.mat <- norm.mat[order(rownames(norm.mat)),]
  ## Converting to a df
  norm.df <- data.frame(Feature = rownames(norm.mat), norm.mat, check.names = FALSE)
  colnames(norm.df)[1] <- feature.out
  ## Writing out df
  if (!is.null(out.dir)) write.table(norm.df, file = paste0(out.dir, '/normalized_counts_', feature.out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  ## Return
  return(norm.mat)
}
## Recount3

## As this function relies on recount3, it requires an internet connection.
r3.dl <- function(species = 'human', project.id = 'BRCA', data.type = 'gene', annotation = 'genecode_v29', out.dir = getwd(), return = FALSE) {
  library(recount3)
  spe.projects <- recount3::available_projects(organism = species)
  ## Checking if the project exists
  if (!project.id %in% spe.projects$project) stop('Provided project.id does not exist !')
  ## Checking if annotation exists
  if (!annotation %in% recount3::annotation_options(species)) stop('Provided annotation does not exist for the provided species !')
  ## Retrieve the data
  my.project <- recount3::create_rse(project_info = spe.projects[spe.projects$project == project.id,], type = data.type, annotation = annotation)
  ## Normalize with VST
  norm.project <- DESeq2::vst(DESeq2::DESeqDataSet(se = my.project, design = as.formula("~1")), blind = TRUE)
  norm.mat <- SummarizedExperiment::assay(norm.project)
  norm.df <- data.frame(Feature = rownames(norm.mat), norm.mat)
  rm(norm.mat)
  ## Saving the RangedSummarizedExperiment object
  if(!is.null(out.dir)) {
    dir.create(path = out.dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(my.project, file = paste0(out.dir, '/', paste(c(species, project.id, data.type, annotation), collapse = '_'), '.RDS'), compress = 'bzip2')
    ## Saving the annotation table
    write.table(SummarizedExperiment::colData(my.project), file = paste0(out.dir, '/', paste(c(species, project.id, 'annotations'), collapse = '_'), '.txt'), sep = "\t", quote = FALSE, row.names = FALSE)
    ## Saving the VST normalized counts table
    write.table(norm.df, file = paste0(out.dir, '/', paste(c(species, project.id, 'vst.normalized.counts'), collapse = '_'), '.txt'), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
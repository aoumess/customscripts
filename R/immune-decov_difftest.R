## PARAMETERS

### B20067_FXDA_01
id.dir <- '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P30_AUMA/P30_AUMA_immune'
annot.file <- '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P30_AUMA/P30_AUMA_immune/P30_AUMA_salmon_immune_annotation.tsv'
annot.colnames <- c('Timepoint', 'C3_response', 'BOR', 'Objective_response', 'Disease_control_benefit', 'Timepoint_OR', 'Timepoint_DCB', 'C2J1_OR', 'Screening_OR', 'C2J1_DCB', 'Screening_DCB')

### B21002_DELE_01
id.dir <- '/home/job/WORKSPACE/B21002_DELE_01/B21002_DELE_01_immune'
annot.file <- '/home/job/WORKSPACE/B21002_DELE_01/B21002_DELE_01_samples_annotation_20210220.csv'
annot.colnames <- c('Macrophage_type', 'F2020', 'Macrophage_type_F2020', 'Basal_F2020', 'Proinf_F2020', 'Treated_Macrophage_type', 'Untreated_Macrophage_type')


## Listing immune-decov result files
id.files <- list.files(path = id.dir, pattern = 'fractions.rds', all.files = TRUE, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
id.tools <- basename(dirname(id.files))
id.list <- setNames(id.files, id.tools)
rm(id.files, id.tools)

## Loading annotations
annot.df <- read.table(annot.file, header = TRUE, sep = "\t")

## Looping !
for (id in names(id.list)) {
  out.dir <- paste0(dirname(id.list[[id]]), '/annot_comparison')
  dir.create(out.dir, recursive = TRUE)
  
  ## Loading deconv results
  decov.res <- readRDS(id.list[[k]])
  decov.mat <- as.matrix(decov.res[,-1])
  rownames(decov.mat) <- decov.res[,1]
  rm(decov.res)
  colnames(decov.mat) <- gsub(pattern = 'X|_R1_trim.fq*', replacement = '', x = colnames(decov.mat))
  ## Sorting according to sample name
  decov.mat <- decov.mat[,order(colnames(decov.mat))]
  ## Checking if our matrix and the annotations are synched
  all(colnames(decov.mat) == annot.df$SampleID)
  
  for (test.term in annot.colnames) {
    ### Performing a KW test (to handle more than 2 classes. Seems very similar to Wilcoxon, with slightly less stringency on raw p-values)
    test.res <- sapply(seq_len(nrow(decov.mat)), function(x) { kruskal.test(decov.mat[x,] ~ annot.df[[test.term]])}, simplify = FALSE)
    ### Converting resultsto a df
    test.df <- data.frame(CellType = rownames(decov.mat), stringsAsFactors = FALSE)
    test.df$Test.score <- vapply(seq_along(test.res), function(x){test.res[[x]]$statistic}, .1)
    test.df$P.value <- vapply(seq_along(test.res), function(x){test.res[[x]]$p.value}, .1)
    ### Adjusting p-values (BH)
    test.df$Adj.p <- p.adjust(test.df$P.value, method = 'BH')
    ### Computing per-class median values
    med.mat <- matrix(nrow = nrow(decov.mat), ncol = nlevels(annot.df[[test.term]]), dimnames = list(rownames(decov.mat), levels(annot.df[[test.term]])))
    for(x in seq_len(nrow(med.mat))) { med.mat[x,] <- aggregate(decov.mat[x,] ~ annot.df[[test.term]], FUN = median)[,2] }
    colnames(med.mat) <- paste0('MEDIAN.', colnames(med.mat))
    ### Computing per-class sd values
    sd.mat <- matrix(nrow = nrow(decov.mat), ncol = nlevels(annot.df[[test.term]]), dimnames = list(rownames(decov.mat), levels(annot.df[[test.term]])))
    for(x in seq_len(nrow(sd.mat))) { sd.mat[x,] <- aggregate(decov.mat[x,] ~ annot.df[[test.term]], FUN = sd)[,2] }
    colnames(sd.mat) <- paste0('SD.', colnames(sd.mat))
    ### Pop size
    popsize <- aggregate(decov.mat[x,] ~ annot.df[[test.term]], FUN = length)[,2]
    test.df <- cbind(test.df, med.mat, sd.mat)
    test.df[paste0('SAMPLES.', levels(annot.df[[test.term]]))] <- rep(popsize, each = nrow(test.df))
    write.table(x = test.df, file = paste0(out.dir, '/', id, '_', test.term, '_KW.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
  }
}


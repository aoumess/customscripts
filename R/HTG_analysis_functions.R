## FUNCTIONS NEEDED FOR THE HTG QC/ANALYSIS WORKFLOW
### XLS data import from multiple files to a matrix
htgxls.import <- function(xls.files = NULL, n.samples = c(24), samplenames.row = c(10), data.startline = c(12)) {
  if (length(xls.files) > 1) {
    if (length(n.samples) == 1) n.samples <- rep(n.samples, length(xls.files))
    if (length(samplenames.row) == 1) samplenames.row <- rep(samplenames.row, length(xls.files))
    if (length(data.startline) == 1) data.startline <- rep(data.startline, length(xls.files))
  }
  ## Import file(s)
  xlz <- sapply(seq_along(xls.files), function(x) { .xl2htg(xls.file = xls.files[x], samplenames.row = samplenames.row[x], n.samples = n.samples[x], data.startline = data.startline[x]) }, simplify = FALSE)
  
  ## Formatting
  # library(dplyr)
  `%>%` <- dplyr::`%>%`
  xlz2 <- xlz %>% purrr::map(tibble::rownames_to_column) %>% purrr::reduce(dplyr::left_join, by = "rowname")
  xlz2 <- as.data.frame(xlz2, row.names = rownames(xlz[[1]]))
  rm(xlz)
  rownames(xlz2) <- xlz2[,1]
  return(as.matrix(xlz2[,-1]))
}

## Internal function to read a single HTG MS Excel file
.xl2htg <- function(xls.file = NULL, n.samples = 24, samplenames.row = 10, data.startline = 12) {
  message('Reading [', xls.file, '] ...')
  rowmax <- suppressMessages(nrow(readxl::read_excel(path = xls.file, sheet = 1, progress = FALSE)))
  last.col <- letters[n.samples + 1]
  xl.title <- suppressMessages(unlist(readxl::read_excel(path = xls.file, sheet = 1, range = paste0('B', samplenames.row, ':', last.col, samplenames.row), col_names = FALSE, progress = FALSE)))
  xl.gnames <- suppressMessages(unlist(readxl::read_excel(path = xls.file, sheet = 1, range = paste0('A', data.startline, ':A', rowmax), col_names = FALSE, progress = FALSE)))
  xl.df <- suppressMessages(readxl::read_excel(path = xls.file, sheet = 1, range = paste0('B', data.startline, ':', last.col, rowmax), col_names = FALSE, progress = FALSE))
  dimnames(xl.df) <- list(xl.gnames, xl.title)
  return(xl.df)
}

### Compute sparsity level
sparsity.level <- function(x = NULL, return.value = TRUE) {
  scprod <- prod(dim(x))
  # message('Total expressed feature count :')
  scZ <- sum(sparseMatrixStats::colCounts(x = x, value = 0))
  # message('\t', scprod - scZ)
  # message('Sparsity level :')
  splev <-  scZ / scprod 
  # message('\t', sprintf('%.5f', splev * 100), '%')
  if(return.value) return(splev)
}

### Get PCA medoids
get.pca.medoid <- function(myvec = NULL, splitvec=NULL, dim = 1) {
  vapply(unique(splitvec), function(x) { median(myvec[splitvec == x]) }, .1)
}

## Convert counts to log10(+1)
int2l10 <- function(x = NULL, epsilon = 1) { log10(x + epsilon)}

## Converts a vector to a color-numeric vector
vec2col <- function(x = NULL) {
  tmp_fac <- as.factor(x)
  tmp_num <- as.numeric(tmp_fac)
  tmp_out <- setNames(object = tmp_num, nm = levels(tmp_fac)[tmp_num])
  return(tmp_out)
  # return(as.numeric(as.factor(x))+1)
}

## Context-colored boxplots
htg.boxplot <- function(matrix.list = NULL, annot.df = NULL, col.item = NULL, col = NULL) {
  oripar <- par()
  par(xaxs = 'i')
  colorder <- order(col)
  # colorder <- seq_along(col)
  annot.df <- annot.df[colorder,]
  legendvec <- sapply(split(x = col,f = names(col), drop = TRUE), unique)
  for(x in seq_along(matrix.list)) {
    matrix.list[[x]] <- matrix.list[[x]][,colorder]
    boxplot(matrix.list[[x]], col = col[colorder], pch = 20, main = paste0(names(matrix.list)[x], ', colored by ', col.item), ylab = 'log10(counts+1)', xaxs = 'i', xaxt = 'n')
    legend(x = ncol(matrix.list[[x]]), y = max(matrix.list[[x]], na.rm = TRUE), legend = names(legendvec), fill = legendvec, xjust = 1)
  }
  suppressWarnings(par(oripar))
}

## Euclidean-forced wrapper for pvclust
dist2.euclidean <- function(x = NULL, ...) {
  amap::Dist(x = t(x), method = 'euclidean', ...)
}

## Pearson-forced wrapper for pvclust
dist2.pearson <- function(x = NULL, ...) {
  amap::Dist(x = t(x), method = 'pearson', ...)
}

## Spearman-forced wrapper for pvclust
dist2.spearman <- function(x = NULL, ...) {
  amap::Dist(x = t(x), method = 'spearman', ...)
}

## Manhattan-forced wrapper for pvclust
dist2.manhattan <- function(x = NULL, ...) {
  amap::Dist(x = t(x), method = 'manhattan', ...)
}

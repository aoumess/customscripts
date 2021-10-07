## FUNCTIONS NEEDED FOR THE HTG QC/ANALYSIS WORKFLOW
### XLS data import function
htgxls.import <- function(xls.file, n.samples = 24, samplenames.row = 10, data.startline = 12) {
  rowmax <- nrow(readxl::read_xls(path = xls.file, sheet = 1))
  last.col <- letters[n.samples + 1]
  xl.title <- unlist(readxl::read_xls(path = xls.file, sheet = 1, range = paste0('B', samplenames.row, ':', last.col, samplenames.row), col_names = FALSE))
  xl.gnames <- unlist(readxl::read_xls(path = xls.file, sheet = 1, range = paste0('A', data.startline, ':A', rowmax), col_names = FALSE ))
  xl.df <- readxl::read_xls(path = xls.file, sheet = 1, range = paste0('B', data.startline, ':', last.col, rowmax), col_names = FALSE)
  dimnames(xl.df) <- list(xl.gnames, xl.title)
  return(xl.df)
}

### Get PCA medoids
get.pca.medoid <- function(myvec = NULL, splitvec=NULL, dim = 1) {
  vapply(unique(splitvec), function(x) { median(myvec[splitvec == x]) }, .1)
}

## Convert counts to log10(+1)
int2l10 <- function(x = NULL, epsilon = 1) { log10(x + epsilon)}

## Converts a vector to a color-numeric vector
vec2col <- function(x = NULL) {
  return(as.numeric(as.factor(x))+1)
}

## Context-colored boxplots
htg.boxplot <- function(matrix.list = NULL, annot.df = NULL, col.item = NULL, col = NULL) {
  oripar <- par()
  par(xaxs = 'i')
  colorder <- order(col)
  annot.df <- annot.df[colorder,]
  for(x in seq_along(matrix.list)) {
    matrix.list[[x]] <- matrix.list[[x]][,colorder]
    boxplot(matrix.list[[x]], border = col[colorder], pch = 20, main = paste0(names(matrix.list)[x], ', colored by ', col.item), ylab = 'log10(counts+1)', xaxs = 'i')
  }
  par(oripar)
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

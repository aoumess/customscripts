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



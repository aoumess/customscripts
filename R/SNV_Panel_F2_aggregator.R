## This script performs the aggregation of "F2" tables (from the Agilent SureSelect XTHS / HaloPlex HS pipelines, corresponding to tabular outputs from Varscan2+Annovar), and additional filtering on variant frequency, ExAC_ALL frequency, and refGene functions.

#####

## Perform the aggregation of F2 output files (from the SSXTHS/HPHS pipeline)
## f2.files           vector(character)       Path(s) to *F2.txt tables
## out.dir            character               Path to output the aggregated data (if NULL, no output written on disk) 
f2.aggreg <- function(f2.files = list.files(path = getwd(), pattern = ".*F2.txt", full.names = TRUE, recursive = TRUE), out.dir = NULL) {
  ## Getting sample names
  sample.names <- unname(vapply(basename(f2.files), function(x) { unlist(strsplit(x = x, split = '_Locatit'))[1]}, 'a'))
  ## Reading header
  f2.headers <- sapply(f2.files, function(x) { read.table(file = x, nrow = 1, header = FALSE, sep = "\t", stringsAsFactors = FALSE) }, simplify = FALSE)
  unique(unlist(Map(length, f2.headers)))
  ## Reading data
  f2.data <- sapply(seq_along(f2.files), function(x) {
    message(paste0('Reading ', sample.names[x], ' ...'))
    f2.rn <- basename(dirname(dirname(f2.files[x])))
    my.df <- try(read.table(file = f2.files[x], skip = 1, header = FALSE, sep = "\t", stringsAsFactors = FALSE, na.strings = '.'))
    if(!is(my.df, class2 = 'try-error')) return(cbind(RunName = f2.rn, Sample = sample.names[x], my.df)) else return(invisible(NULL))
  }, simplify = FALSE)
  ## Removing empty samples
  empty.samples <- sapply(f2.data, is.null)
  f2.headers <- f2.headers[!empty.samples]
  f2.data <- f2.data[!empty.samples]
  f2.merged <- Reduce(rbind, f2.data)
  ## Correcting bad last cols
  coldiff <- ncol(f2.merged) -unique(unlist(Map(length, f2.headers)))
  f2.form <- f2.merged[,-c((ncol(f2.merged)-(coldiff*2)+3):(ncol(f2.merged)-(coldiff*1)))]
  ## Adding colnames
  colnames(f2.form) <- c('RunName', 'Sample', unlist(f2.headers[[1]], use.names = FALSE))
  ## Fixing freq format
  f2.form$FREQ <- as.numeric(sub(pattern = '%', replacement = '', x = f2.form$FREQ)) / 100
  ## Dumping results (when needed)
  if (!is.null(out.dir)) write.table(f2.form, file = paste0(out.dir, '/F2_merged.txt'), sep = "\t", quote = FALSE, row.names = FALSE)
  ## Func output
  return(f2.form)
}


## Perform the filtering of the f2.aggreg output
## aggreg.df          data.frame          df output from f2.aggreg()
## freq.min           numeric             minimal variant freq to retain
## ExALL.max          numeric             maximal variant freq in ExAC_ALL to retain
## refGen.funcs       vector(character)   refGen functions to retain
## out.dir            character       Path to output the aggregated data (if NULL, no output written on disk)
f2.filter <- function(aggreg.df = NULL, freq.min = .01, ExALL.max = .005, refGen.funcs = c('exonic', 'splicing'), out.dir = NULL) {
  ## Cut 
  filt.df <- aggreg.df[which(((aggreg.df$ExAC_ALL < ExALL.max) | is.na(aggreg.df$ExAC_ALL)) & aggreg.df$FREQ > freq.min & aggreg.df$Func.refGene %in% refGen.funcs), ]
  if (!is.null(out.dir)) write.table(filt.df, file = paste0(out.dir, '/F2_filtered_Freq', freq.min, '_ExALL', ExALL.max, '_', paste(refGen.funcs, collapse = '.'), '.txt'), sep = "\t", quote = FALSE, row.names = FALSE)
  return(filt.df)
}


#####

## USAGE
# f2.dir <- '/home/job/WORKSPACE/B19033_JEMI_01/ANALYSIS/P31_JBMI_HaloplexHS'
# f2.files <- list.files(path = f2.dir, pattern = '*F2.txt$', full.names = TRUE, recursive = TRUE)
# foo <- f2.aggreg(f2.files = f2.files, out.dir = f2.dir)
# bar <- f2.filter(aggreg.df = foo, out.dir = f2.dir)

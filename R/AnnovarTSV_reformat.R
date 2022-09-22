## This script performs the aggregation then filtering of Annovar-annotated TSV files  tables (from the Agilent SureSelect XTHS / HaloPlex HS pipelines, corresponding to tabular outputs from Varscan2/FreeBayes+Annovar), and additional filtering on variant frequency, ExAC_ALL frequency, and refGene functions.

annovar.tsv_aggreg <- function(annovar.tsv.files = list.files(path = getwd(), pattern = ".*_multianno.txt", full.names = TRUE, recursive = TRUE), out.dir = NULL, return.data = FALSE) {
  ## Getting sample names
  sample.names <- unname(vapply(basename(annovar.tsv.files), function(x) { unlist(strsplit(x = x, split = '_Locatit'))[1]}, 'a'))
  ## Reading header
  at.headers <- sapply(annovar.tsv.files, function(x) { read.table(file = x, nrow = 1, header = FALSE, sep = "\t", stringsAsFactors = FALSE) }, simplify = FALSE)
  ## Checking they all have the same nomber of columns
  unique(unlist(Map(length, at.headers)))
  ## Reading data
  at.data <- sapply(seq_along(annovar.tsv.files), function(x) {
    message(paste0('Reading ', sample.names[x], ' ...'))
    at.rn <- basename(dirname(annovar.tsv.files[x]))
    my.df <- try(read.table(file = annovar.tsv.files[x], skip = 1, header = FALSE, sep = "\t", stringsAsFactors = FALSE, na.strings = '.'))
    if(!is(my.df, class2 = 'try-error')) return(cbind(RunName = at.rn, Sample = sample.names[x], my.df)) else return(invisible(NULL))
  }, simplify = FALSE)
  ## Removing empty samples
  empty.samples <- sapply(at.data, is.null)
  at.headers <- at.headers[!empty.samples]
  at.data <- at.data[!empty.samples]
  at.merged <- Reduce(rbind, at.data)
  rm(at.data)
  ## Remove bad columns
  at.merged <- at.merged[,-c((unique(unlist(Map(length, at.headers)))+2):(ncol(at.merged)-2))]
  colnames(at.merged) <- c('RunName', 'Sample', unlist(at.headers[[1]], use.names = FALSE), 'VALZ')
  rm(at.headers)
  # write.table(x = at.merged, file = '/home/job/WORKSPACE/PROJECTS/B21031_JEMI_01/TEST_RFILTER/at.merged.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
  ## Split INFO fields
  iname.head <- unlist(strsplit(x = at.merged[1, (ncol(at.merged)-1)], split = ':'))
  ivalz.split <- sapply(seq_len(nrow(at.merged)), function(l) unlist(strsplit(x = at.merged[l, ncol(at.merged)], split = ':')))
  rownames(ivalz.split) <- iname.head
  at.form <- cbind(at.merged[, -c((ncol(at.merged)-1):ncol(at.merged))], t(ivalz.split))
  rm(at.merged)
  # write.table(x = at.form, file = '/home/job/WORKSPACE/PROJECTS/B21031_JEMI_01/TEST_RFILTER/at.form.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
  ## Fixing AD
  at.form$AD <- vapply(at.form$AD, function(x) as.integer(rev(unlist(strsplit(x = x, split = ',')))[1]), 1L)
  ## Converting some info to numeric
  at.form$DP <- as.integer(at.form$DP)
  ## Computing FREQ
  at.form$FREQ <- at.form$AD / at.form$DP
  ## Dump
  if(!is.null(out.dir)) write.table(x = at.form, file = paste0(out.dir, '/ALL_annotated_variants_unfiltered.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  ## Return
  if (return.data) return(at.form)
}

## Perform the filtering of the f2.aggreg output
## aggreg.df          data.frame          df output from annovar.tsv_aggreg()
## filter.list        list                list of terms to evaluate for filtering (see example)
## out.dir            character|NULL      Path to output the aggregated data (if NULL, no output written on disk)
## return.data        bool                Function returns filtered df
annovar.tsv_filter <- function(aggreg.df = NULL, filter.list = NULL, out.dir = getwd(), return.data = FALSE) {
  
  if (!is.null(filter.list)) {
    message(nrow(aggreg.df))
    message('FILTERING ...')
    for (x in names(filter.list)) {
      message(paste0('\t', x))
      if ('keep' %in% names(filter.list[[x]])) aggreg.df <- if (filter.list[[x]]$keepNA) aggreg.df[is.na(aggreg.df[[x]]) | grepl(pattern = paste(filter.list[[x]][['keep']], collapse = '|'), x = aggreg.df[[x]]), ] else aggreg.df[grep(pattern = paste(filter.list[[x]][['keep']], collapse = '|'), x = aggreg.df[[x]]), ]
      # dim(tsv_df)
      if ('discard' %in% names(filter.list[[x]])) aggreg.df <- if (filter.list[[x]]$keepNA) aggreg.df[is.na(aggreg.df[[x]]) | !grepl(pattern = paste(filter.list[[x]][['discard']], collapse = '|'), x = aggreg.df[[x]]), ] else aggreg.df[grep(pattern = paste(filter.list[[x]][['keep']], collapse = '|'), x = aggreg.df[[x]], invert = TRUE), ]
      # dim(tsv_df)
      if ('eval' %in% names(filter.list[[x]])) aggreg.df <- if (filter.list[[x]]$keepNA) aggreg.df[is.na(aggreg.df[[x]]) | eval(parse(text=paste0('aggreg.df[[x]] ', filter.list[[x]][[1]]))), ] else aggreg.df[!is.na(aggreg.df[[x]]) & eval(parse(text=paste0('aggreg.df[[x]] ', filter.list[[x]][[1]]))), ]
      message(nrow(aggreg.df))
    }
  
    ## Converting filter.list to text
    filt.txt <- paste(vapply(names(filter.list), function(x) paste(c(x, names(filter.list[[x]]), unlist(filter.list[[x]])), collapse = '.'), 'a'), collapse = '+')
    filt.txt <- gsub(pattern = ' ', replacement = '_', x = filt.txt)
    aggreg.df$FILTER <- filt.txt
    ## Dumping
    if (!is.null(out.dir)) write.table(aggreg.df, file = paste0(out.dir, '/ALL_annotated_variants_filtered.tsv'), sep = "\t", quote = FALSE, row.names = FALSE)
    ## Returning
    if(return.data) return(aggreg.df)
  }
}


## VARZ
workdir = '/home/job/WORKSPACE/PROJECTS/B21031_JEMI_01/TEST_RFILTER'
filter.list <- list(
  ExAC_ALL = list(eval = '<0.005', keepNA = TRUE)
  , Func.refGene = list(keep = c('exonic', 'splicing'), keepNA = TRUE)
  , ExonicFunc.refGene = list(discard = c('synonymous SNV'), keepNA = TRUE)
  # , FREQ = list(eval = '<.001', keepNA = TRUE)
)

## RUN
setwd(workdir)
at.list <- list.files(path = getwd(), pattern = ".*_multianno.txt", full.names = TRUE, recursive = TRUE)
## Aggregating / formatting
at.all <- annovar.tsv_aggreg(annovar.tsv.files = at.list, out.dir = workdir, return.data = TRUE)
## Filtering
at.filt <- annovar.tsv_filter(aggreg.df = at.all, filter.list = filter.list, out.dir = workdir, return.data = TRUE)

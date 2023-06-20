### Retrieve genomic annotation from UCSC (using a UCSC genome name) and annotate a regional dataframe

### REQUIREMENTS :
# CRAN : dplyr, tibble
# BioConductor : GenomeInfoDb, GenomicRanges
# Github : jeffbhasin/goldmine

## Example call on a DF with columns chr/start/end :
## 1) get annotation for a defined genome and store them in cache
##  my.annot.obj <- get.annotation(genome = 'hg19', cache.dir = '.')
## 2) annotate a genomic region DF
##  my.annot.df <- regions.annotate(regions.df = raw.df, annot.obj = my.annot.obj)

### Get genomic annotations ====
get.annotation <- function(genome = NULL, annot.type = c('cytoBandIdeo', 'cpgIslandExt', 'genomicSuperDups', 'dgvMerged'), cache.dir = NULL, sync = FALSE, save.blob = FALSE) {
  
  ## Getting the list of available registered genomes @ UCSC (from GenomeInfoDb)
  avail.genz <- sort(GenomeInfoDb::registered_UCSC_genomes()$genome)
  
  if (is.null(genome)) stop("A UCSC genome name is required. It should be one of : '", paste(avail.genz, collapse = "', '"), "'")
  if(is.null(cache.dir)) stop('A cache directory to write files is required !')
  if(!file.exists(cache.dir)) stop(paste0('Cache directory [', cache.dir, '] does not exist !'))
  
  message('An internet connexion is required when cached data are not available.')
  
  genome <- genome[genome %in% avail.genz]
  if(length(genome) == 0) stop('None of the specified genome(s) is available @ UCSC !')
  
  annot.obj <- list()
  message(paste0('Getting annotations for ', genome, ' ...'))
  
  ## Trying to get chromosomes info
  message('Getting chromosomes info ...')
  trycI <- try(cI <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome, assembled.molecules.only = TRUE, as.Seqinfo = TRUE), silent = TRUE)
  ## Getting RefSeq genes
  message('Getting RefSeq genes ...')
  annot.obj$refGene <- goldmine::getGenes(geneset = 'refseq', genome = genome, cachedir = cache.dir, sync = sync)
  
  ## Cycling in annotation tables
  for (tt in annot.type) {
    message('Getting ', tt, ' data ...')
    tryAn <- try(annot.obj[[tt]] <- goldmine::getUCSCTable(table = tt, genome = genome, cachedir = cache.dir, sync = sync), silent = TRUE)
    if (!class(tryAn)[1] == 'try-error') annot.obj[[tt]][['bin']] <- NULL
  }
  ## Convert to GRanges
  message('Converting to GRanges ...')
  for (tt in seq_along(annot.obj)) {
    annot.obj[[tt]] <- GenomicRanges::makeGRangesFromDataFrame(df = annot.obj[[tt]], keep.extra.columns = TRUE, ignore.strand = FALSE, seqnames.field = colnames(annot.obj[[tt]])[1], start.field = colnames(annot.obj[[tt]])[2], end.field = colnames(annot.obj[[tt]])[3], starts.in.df.are.0based = TRUE)
    if (!class(trycI)[1] == 'try-error') {
      GenomeInfoDb::seqlevels(annot.obj[[tt]], pruning.mode = 'coarse') <- GenomeInfoDb::seqlevels(cI)
      GenomeInfoDb::seqinfo(annot.obj[[tt]]) <- cI
    }
  }
  
  ## Merging refGenes
  message('Merging refGenes ...')
  `%>%` <- dplyr::`%>%`
  annot.obj$refGene <- suppressMessages(annot.obj$refGene %>% tibble::as_tibble() %>% dplyr::group_by(seqnames, strand, name) %>% dplyr::summarize(start = min(start), end = max(end)) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE))
  
  ## Saving RDA
  if(save.blob) {
    message('Saving data blob ...')
    save(annot.obj, file = paste0(cache.dir, '/', genome, '_annot.RDA'))
  }
  
  ## Return
  return(annot.obj)
}

### Annotate regions ====
regions.annotate <- function(regions.df = NULL, chr.colname = 'Chr', start.colname = 'Start', end.colname = 'End', annot.obj = NULL, add.width = TRUE, skip.genes = FALSE) {
  
  ## Add width
  # message('COLNAMES : ', paste(colnames(regions.df), collapse = ', '))
  # message('EXPECTED : ', paste(c(chr.colname, start.colname, end.colname), collapse = ', '))
  pos.colnames <- match(c(chr.colname, start.colname, end.colname), colnames(regions.df))
  # message(pos.colnames)
  if(add.width) regions.df <- cbind(regions.df[,1:max(pos.colnames)], Width = regions.df[,pos.colnames[3]] - regions.df[,pos.colnames[2]] + 1, regions.df[,-c(1:max(pos.colnames))])
  
  ## Cytobands
  if (!is.null(annot.obj$cytoBandIdeo)) {
    message('Annotating regions with cytologic bands ...')
    regions.df$Cyto.End <- regions.df$Cyto.Start <- NA
    ## Start
    fin.cstart.gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions.df, seqnames.field = chr.colname, start.field = start.colname, end.field = start.colname, keep.extra.columns = FALSE)
    cyto.start.ol <- GenomicRanges::findOverlaps(query = fin.cstart.gr, subject = annot.obj$cytoBandIdeo)
    regions.df$Cyto.Start[S4Vectors::from(cyto.start.ol)] <- annot.obj$cytoBandIdeo$name[S4Vectors::to(cyto.start.ol)]
    rm(fin.cstart.gr, cyto.start.ol)
    ## End
    fin.cend.gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions.df, seqnames.field = chr.colname, start.field = end.colname, end.field = end.colname, keep.extra.columns = FALSE)
    cyto.end.ol <- GenomicRanges::findOverlaps(query = fin.cend.gr, subject = annot.obj$cytoBandIdeo)
    regions.df$Cyto.End[S4Vectors::from(cyto.end.ol)] <- annot.obj$cytoBandIdeo$name[S4Vectors::to(cyto.end.ol)]
    rm(fin.cend.gr, cyto.end.ol)
  }
  
  ## DGV
  if (!is.null(annot.obj$dgvMerged)) {
    message('Annotating regions with DGV CNV ...')
    fin.gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions.df, seqnames.field = chr.colname, start.field = start.colname, end.field = end.colname, keep.extra.columns = TRUE)
    regions.df$DGV.CNV <- GenomicRanges::countOverlaps(query = fin.gr, subject = annot.obj$dgvMerged)
    rm(fin.gr)
  }
  
  ## GenomicSuperDups
  if (!is.null(annot.obj$genomicSuperDups)) {
    message('Annotating regions with UCSC genomic superduplications ...')
    fin.gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions.df, seqnames.field = chr.colname, start.field = start.colname, end.field = end.colname, keep.extra.columns = TRUE)
    regions.df$GenomicSuperDups <- GenomicRanges::countOverlaps(query = fin.gr, subject = annot.obj$genomicSuperDups)
    rm(fin.gr)
  }
  
  ## CpGislands
  if (!is.null(annot.obj$cpgIslandExt)) {
    message('Annotating regions with UCSC CpG islands ...')
    fin.gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions.df, seqnames.field = chr.colname, start.field = start.colname, end.field = end.colname, keep.extra.columns = TRUE)
    regions.df$CpGislands <- GenomicRanges::countOverlaps(query = fin.gr, subject = annot.obj$cpgIslandExt)
    rm(fin.gr)
  }
  
  ## Refseq Genes
  if (!is.null(annot.obj$refGene) & !skip.genes) {
    message('Annotating regions with RefSeq genes ...')
    fin.gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions.df, seqnames.field = chr.colname, start.field = start.colname, end.field = end.colname, keep.extra.columns = TRUE)
    regions.df$RefSeq.Genes <- GenomicRanges::countOverlaps(query = fin.gr, subject = annot.obj$refGene)
    genes.oldf <- as.data.frame(GenomicRanges::findOverlaps(query = fin.gr, subject = annot.obj$refGene))
    regions.df$Symbols <- NA
    unique.hits <- unique(genes.oldf$queryHits)
    genes.list <- sapply(unique.hits, function(x) {
      paste(annot.obj$refGene$name[genes.oldf$subjectHits[genes.oldf$queryHits == x]], collapse = ',')
    }, simplify = TRUE)
    regions.df$Symbols[unique.hits] <- genes.list
    rm(fin.gr, genes.oldf, unique.hits, genes.list)
  }
  return(regions.df)
}



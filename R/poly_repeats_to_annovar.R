## This scripts gets position of polymer repeats in a BSgenome and outputs these positions as an ANNOVAR DB file.

# setwd("/home/job/WORKSPACE/SSR")


poly_repeats_to_annovar_db <- function(genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", min.size = 5, k.vec = 1:24) {
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)
  knames <- seqnames(BSg.obj)[k.vec]
  
  require(foreach)
  
  poly.repeats <- foreach(k = knames, .combine = "rbind") %do% {
    message(k)
    kseq <- as.character(BSgenome::getSeq(BSg.obj, k))
    ksplit <- unlist(strsplit(kseq, ""))
    rm(kseq)
    kdf <- data.frame(chr=k, unclass(rle(ksplit)))
    rm(ksplit)
    kdf$end <- cumsum(kdf$lengths)
    kdf$start <- c(1L, kdf$end[-nrow(kdf)]+1L)
    kdf$end <- kdf$start
    kdf <- kdf[kdf$lengths >= min.size & kdf$values != "N",]
    kdf$alt <- "-"
    kdf <- kdf[,c(1,5,4,3,6,2)]
    kdf$lengths <- paste0("poly.", kdf$values, ":", kdf$lengths)
    levels(kdf[,1]) <- sub(pattern = "^chr", replacement = "", x = levels(kdf[,1]), ignore.case = TRUE)
    colnames(kdf) <- c("#Chr", "Start", "End", "Ref", "Alt", paste0("Poly.repeats_", min.size))
    return(kdf)
  }
  write.table(poly.repeats, file = paste0(genome, "_poly.repeats_", min.size, "_", format(Sys.Date(), "%Y%m%d"), ".txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  rm(poly.repeats)
}
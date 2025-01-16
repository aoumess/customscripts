
## Generate TCN-CBS from ASCAT results
# setwd("/mnt/data_colibri/share/analyses/Projets/AUAU/P25_AUAU2/ASCAT/ASCAT_pHe0.50_pO0.03_sL5_pen50_nrf0.5_SER/")
# optimdirs <- list.files(path = getwd(), full.names = TRUE, recursive = TRUE, pattern = "_optimal$", include.dirs = TRUE, all.files = FALSE, ignore.case = FALSE)
# foreach(optimgo = optimdirs, .combine = "rbind") %do% {
#   as.res <- readRDS(list.files(path = optimgo, pattern = ".RDS$", full.names = TRUE, recursive = FALSE, include.dirs = FALSE, ignore.case = TRUE))
#   outfile <- paste0(optimgo, "/", as.res$segments$sample[1], "_", basename(optimgo), "_TCN.cbs")
#   outdf <- as.res$segments
#   outdf$chr <- unlist(cs$chrom2chr[paste0("chr", outdf$chr)])
#   outdf$Width <- outdf$endpos - outdf$startpos + 1
#   outdf$TCN <- outdf$nMajor + outdf$nMinor
#   outdf$Ploidy <- as.res$ploidy
#   outdf$GoF <- as.res$goodnessOfFit
#   outdf <- outdf[,c(1:4,7,8,5,6,9,10)]
#   colnames(outdf)[1:4] <- c(as.res$segments$sample[1], "Chr", "Start", "End")
#   maelstrom::write.table.fast(x = outdf, file = outfile)
# }

## Parsing args
require(optparse)
option_list <- list(
  make_option(c("-i", "--in"), type = "character", default= NULL, help="Path leading to a TCN-CBS file"),
  make_option(c("--genome"), type = "character", default="hg19", help="A genome build name, as supported by R package 'chromosomes'")
)
opt <- parse_args(OptionParser(option_list=option_list))
# print (opt)

infile <- opt[["in"]]
genome <- opt[["genome"]]

## Checking args
if (is.null(infile)) stop("A TCN-CBS input file is required !")
if (!file.exists(infile)) stop(paste0("Input file '", infile, "' could not be found !"))

## Loading chromosomes data
data(list = genome, package = "chromosomes", envir = environment())
# data(list = "hg19", package = "chromosomes", envir = environment())

## Loading TCN-CBS file
as.res <- maelstrom::read.table.fast(file = infile)

## Computing scores
ploidy <- as.res$Ploidy[1]
gof <- as.res$GoF[1]
### Sum of segments lengths multiplied by the absolute difference to diploidy, divided by the genome length.
SKOR1 <- sum(abs(as.res$TCN - 2) * (as.res$End - as.res$Start + 1)) / cs$genome.length
### Sum of segments lengths multiplied by the absolute difference to exact ASCAT2-predicted ploidy, divided by the genome length.
SKOR2 <- sum(abs(as.res$TCN - ploidy) * (as.res$End - as.res$Start + 1)) / cs$genome.length
### Sum of segments lengths multiplied by the absolute difference to rounded ASCAT2-predicted ploidy, divided by the genome length.
SKOR3 <- sum(abs(as.res$TCN - round(ploidy)) * (as.res$End - as.res$Start + 1)) / cs$genome.length
### Nb of segments
NBSEG <- nrow(as.res)
### Compute a score on a per-chromosome basis, taking the longest ACN (in bp) as normal basis
require(foreach)
LOKAL1 <- sum(foreach (k = unique(as.res$Chr), .combine = "c") %do% {
  require(dplyr)
  miniseg <- as.tbl(as.res[as.res$Chr == k,])
  minisegtbl <- as.tbl(miniseg)
  basistbl <- as.data.frame(minisegtbl %>% group_by(TCN) %>% summarise(sum(Width)))
  basis <- basistbl[which.max(basistbl[,2]),1]
  KSCOR <- sum(abs(miniseg$TCN  - basis) * miniseg$Width) / cs$genome.length
  return(KSCOR)
})

## Generating output
outdf <- data.frame(SCORE1 = SKOR1, SCORE2 = SKOR2, SCORE3 = SKOR3, NBSEG = NBSEG, LOCAL1 = LOKAL1, stringsAsFactors = FALSE)
outfile <- paste0(sub(pattern = "\\.cbs$", replacement = "_SCORED.txt", x = infile, ignore.case = TRUE))
maelstrom::write.table.fast(x = outdf, file = outfile)

## END

# pdf("SCORES.pdf", width = 29.7/cm(1), height = 21/cm(1))
# plot(SKORZ$SCORE1, pch = 20, col = "black", type = "b", main = "Scores\n(B:v1 / R:v2 / G:v3 / P:Local / O:min-log10(NbSeg))", xlab = "Samples", ylab = "Score")
# lines(SKORZ$SCORE2, pch = 20, col = "red", type = "b")
# lines(SKORZ$SCORE3, pch = 20, col = "green", type = "b")
# lines(SKORZ$LOCAL1, pch = 20, col = "purple", type = "b")
# lines(log10(SKORZ$NbSeg)-min(log10(SKORZ$NbSeg)), pch = 20, col = "orange", type = "b")
# abline(h = 0, lty = 3)
#
# plot(density(SKORZ$LOCAL1, bw = .1), col = "purple", main = "Score densities\n(B:v1 / R:v2 / G:v3 / P:Local / O:min-log10(NbSeg))")
# lines(density(SKORZ$SCORE1, bw = .1), col = "black")
# lines(density(SKORZ$SCORE2, bw = .1), col = "red")
# lines(density(SKORZ$SCORE3, bw = .1), col = "green")
# lines(density(log10(SKORZ$NbSeg)-min(log10(SKORZ$NbSeg)), bw = .1), col = "orange")
#
# dev.off()



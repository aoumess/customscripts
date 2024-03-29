## This script generates a table containing GIS (genomic instability scores) from an EaCoN TCN output using the ASCAT segmenter.

## Requirements :
## - CRAN : optparse, devtools, foreach, dplyr
## - GITHUB : https://github.com/aoumess/chromosomes

## Usage :
## R CMD BATCH EaCoN_TCN_GIS_autoscorer -i path/to/eacon_results/samplename/ -g hg19 -s ASCAT

## Parsing args
option_list <- list(
  optparse::make_option(c("-i", "--input-dir"), type = "character", default= NULL, help="Path leading to an EaCoN result directory (should end with a samplename : /path/to/samplename/)"),
  optparse::make_option(c("-g", "--genome"), type = "character", default="hg19", help="A genome build name, as supported by R package 'chromosomes'"),
  optparse::make_option(c("-s", "--segmenter"), type = "character", default="ASCAT", help="A segmenter name, supported by EaCoN (should be 'ASCAT', 'SEQUENZA' or 'FACETS')")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
print (opt)

indir <- opt[["input-dir"]]
genome <- opt[["genome"]]
segmenter <- opt[["segmenter"]]

## Checking args
if (is.null(indir)) stop("An EaCoN results directory is required !")
if (!file.exists(indir)) stop(paste0("Input directory '", infile, "' could not be found !"))

## Checking if chromosomes is installed
if (!('chromosomes' %in% installed.packages())) {
  if (!('devtools' %in% installed.packages())) {
    install.packages('devtools')
  }
  devtools::install_github('aoumess/chromosomes')
}

## Loading chromosomes data
data(list = genome, package = "chromosomes", envir = environment())

## Loading the gammaEval file
samplename <- basename(indir)
tcnsegdir <- paste0(indir, '/', segmenter, '/ASCN/')
gE <- read.table(paste0(tcnsegdir, samplename, '.gammaEval.txt'), header = TRUE, sep = "\t")

## Looking for the optimal model
gE$round.psi <- round(gE$psi)
gE$round.psi.diff <- abs(gE$psi - gE$round.psi)
gE$gof.rank <- rank(gE$GoF)
gE$rpd.rank <- rank(gE$round.psi.diff)
gE$score <- gE$gof.rank - gE$rpd.rank

best.gamma <- sprintf("%.2f", gE$gamma[which.max(gE$score)])

## Loading TCN-CBS file
as.res <- readRDS(file = paste0(tcnsegdir, '/gamma', best.gamma, '/', samplename, '.ASCN.ASCAT.RDS'))

## Computing needed values
as.res$segments$TCN <- as.res$segments$nMajor + as.res$segments$nMajor
as.res$segments$Width <- as.res$segments$endpos - as.res$segments$startpos + 1

## Computing scores
ploidy <- as.res$ploidy$ascat
gof <- as.res$goodnessOfFit
### Sum of segments lengths multiplied by the absolute difference to diploidy, divided by the genome length.
SKOR1 <- sum(abs(as.res$segments$TCN - 2) * as.res$segments$Width) / cs$genome.length
### Sum of segments lengths multiplied by the absolute difference to exact ASCAT2-predicted ploidy, divided by the genome length.
SKOR2 <- sum(abs(as.res$segments$TCN - ploidy) * as.res$segments$Width) / cs$genome.length
### Sum of segments lengths multiplied by the absolute difference to rounded ASCAT2-predicted ploidy, divided by the genome length.
SKOR3 <- sum(abs(as.res$segments$TCN - round(ploidy)) * as.res$segments$Width) / cs$genome.length
### Nb of segments
NBSEG <- nrow(as.res$segments)
### Compute a score on a per-chromosome basis, taking the longest ACN (in bp) as normal basis
library(foreach)
suppressPackageStartupMessages(library(dplyr))
  LOKAL1 <- sum(foreach (k = unique(as.res$segments$chr), .combine = "c") %do% {
    miniseg <- tibble::as_tibble(as.res$segments[as.res$segments$chr == k,])
    basistbl <- miniseg %>% group_by(TCN) %>% summarise(sum(Width))
    basis <- basistbl$TCN[which.max(as.data.frame(basistbl)[,2])]
    KSCOR <- sum(abs(miniseg$TCN  - basis) * miniseg$Width) / cs$genome.length
    return(KSCOR)
  })

## Generating output
outdf <- data.frame(Samplename = samplename, Best_gamma = best.gamma, Ploidy = ploidy, Cellularity = as.res$aberrantcellfraction, GoF = gof, Psi = as.res$psi, SCORE1 = SKOR1, SCORE2 = SKOR2, SCORE3 = SKOR3, NBSEG = NBSEG, LOCAL1 = LOKAL1, stringsAsFactors = FALSE)
write.table(x = outdf, file = paste0(indir, '/', samplename, '_GIS.txt'), quote = FALSE, sep = "\t", row.names = FALSE)


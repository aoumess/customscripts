## skmeans_run
##
## DESCRIPTION : A wrapper to ease the use of the skmeans package for R. This allows
##               spherical-kmeans clustering using various methods.
##
## AUTHOR:
## Bastien JOB (bastien.job@gustaveroussy.fr)
##
## DEPENDS ON:
## . skmeans, foreach, kmndirs, cluster
## . (gmeans, CLUTO : external softwares - optional)
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 3.3.1
##
## VERSION NOTES:
##
## v2.0 20241008
## . Rewrote the whole skmeans_run() function.
## . Added support for svg+png.
## . Simplified the control of parameters by forwarding the "control" list to skmeans::skmeans.
##
## v1.1 20160426
## . Limited the methods for comparison to six methods (others made comparison fail)
## . Changed default maxIter to 5000 (from 3000)
## . Added a help function
## . Corrected a bug where features scores and contributions to metagenes were not saved

## v1.0 20160302
## . First release


## Just gives the names of built-in skmeans methods
skmeans.methods.list <- function(echo = FALSE) {
  skmlist <- c("genetic", "pclust", "kmndirs", "gmeans", "cluto")
  if (echo) { print('<<NOTA : "kmndirs" requires the kmndirs package installed ; "gmeans" and "cluto" require the corresponding eponymic external softwares installed.>>'); print(skmlist) }
  return(skmlist)
}

## Help for main function
help.skmeans_run <- function() {
  cat("
skmeans_run <- function(data = NULL, k.test = 2:15, method = \"standard\", control = NULL, odir = NULL)

data            :  (matrix)             A features (rows) by samples (cols) matrix of numerical values.
k.test          :  (integer)            A vector of integers corresponding to desired clusters.
method          :  (character)          The skmeans method to use. Run skmeans.methods.list() to obtain a list of available methods.
control         :  (list)               Controls the parameters for each method. See ?skmeans::skmeans.
odir            :  (character)          The directory where results will be put (default NULL keeps current working directory).

")
}
skmeans_run <- function(data = NULL, k.test = 2:5, method = 'standard', control = NULL, odir = getwd()) {
  
  ## Checks
  ### Dirs
  if (is.null(odir)) odir <- getwd()
  if (!dir.exists(odir)) stop(paste0("Can't find output directory : ", odir, " !"))
  ### Methods
  if (!(method %in% skmeans.methods.list(echo = FALSE))) stop("Unsupported method ! Please use a method listed in skmeans.methods.list()")
  #### genetic
  if (tolower(method) == 'genetic') message('Method "gmeans" requires a control list like this : control = list(maxiter = 12, popsize = 6, mutations = .1) where "maxiter" is the max amount of iterations, "popsize" is the population size for the genetic algorithm, and "mutations" the mutation rate at each generation')
  #### pclust
  if (tolower(method) == 'pclust') message('Method "pclust" requires a control list like this : control = list(maxiter = 100, nruns = 1, maxchains = 0) where "maxiter" is the max amount of iterations, "nruns" is the number of fixed-point runs, and "maxchains" the maximal length of the Kernighan-Lin chains')
  #### kmndirs
  if (tolower(method) == 'kmndirs') message('Method "kmndirs" requires the R package "kmndirs" installed, and a control list like this : control = list(maxiter = 10, nstart = 100) where "maxiter" is the max amount of iterations, and "nstart" the number of starting points to compute the starting value for the iteration stage')
  #### gmeans
  if (tolower(method) == 'gmeans') message('Method "gmeans" requires the gmeans binary and a control list like this : control = list(gmeans = "/path/to/gmeans", control = paste0(" -a ", gmeans.a)), where "gmeans" is the path to the gmeans binary, "control" is a character chain giving additional parameters, and "gmeans.a" the gmeans submethod [s|e|b|k|d]')
  #### CLUTO
  if (tolower(method) == "cluto") message('Method "cluto" requires the CLUTO vcluster binary and a control list like this : control = list(vcluster = "/path/to/vcluster", colmodel = "", control = paste0(" -clmethod=", clutomethod)), where "vcluster" is the path to the CLUTO vcluster binary, "control" is a character chain giving additional parameters, and "clutomethod" the cluto submethod [rbr|rb|direct|agglo|graph|bagglo]')
  
  ## Define methodword
  if (tolower(method %in% 'gmeans')) methodword <- paste0(method, '-', sub(pattern = ' -a ', replacement = '', x = control[['control']])) else if (tolower(method %in% 'cluto')) methodword <- paste0(method, '-', sub(pattern = ' -clmethod=', replacement = '', x = control[['control']])) else methodword <- method
  
  ## Run skmeans
  skmeans.res <- lapply(k.test, function(k) {
    message(paste0("Testing k=", k, " ..."))
    skmeans::skmeans(x = t(data), method = method, k = k, control = control)
  })
  names(skmeans.res) <- k.test
  
  ## Create output dir
  outdir <- paste0(odir, '/', methodword)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  ## Output rootname
  or <- paste0(outdir, '/', paste(c('skmeans', methodword, paste(c('k', range(k.test, na.rm = TRUE)), collapse = '.'), paste0('mI.', maxIter)), collapse = '_'))
  
  ## Save result object
  saveRDS(skmeans.res, paste0(or, '_results.RDS'))
  
  ## Plotting silhouettes
  silh.mean <- silh.q25 <- silh.q75 <- vector()
  pdf(file = paste0(or, '_silhouettes.pdf'), width = 29.7/cm(1), height = 21/cm(1))
  for (k in 1:length(skmeans.res)) {
    mysil <- cluster::silhouette(skmeans.res[[k]])
    silh.mean <- c(silh.mean, mean(mysil[,3], na.rm = TRUE))
    silh.q25 <- c(silh.q25, quantile(mysil[,3], .25, na.rm = TRUE))
    silh.q75 <- c(silh.q75, quantile(mysil[,3], .75, na.rm = TRUE))
    plot(mysil, main = paste0("Silhouettes for k = ", names(skmeans.res[k])))
    rm(mysil)
  }
  dev.off()
  ## Silhouettes summary
  svg(filename = paste0(or, '_silhouettes_summary.svg'), width = 1280/96, height = 1024/96)
  plot(as.numeric(names(skmeans.res)), as.numeric(silh.q75), type = "b", pch = 20, main = "All silhouettes", xlab = "k", ylab = "Silhouette (Q25, mean, Q75)", col = 2, ylim = range(c(silh.q75, silh.q25, silh.mean), finite = TRUE, na.rm = TRUE))
  lines(as.numeric(names(skmeans.res)), silh.mean, type = "b", pch = 20)
  lines(as.numeric(names(skmeans.res)), silh.q25, type = "b", pch = 20, col = 2)
  abline(h = 0, lty = 2)
  svg_off()
  ## Save silhouettes data
  kout_df <- data.frame(k = k.test, silhouette = silh.mean, stringsAsFactors = FALSE)
  write.table(x = kout_df, file = paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_silhouettes.txt"), sep="\t", quote = FALSE, row.names = FALSE)
  WriteXLS::WriteXLS(x = kout_df, ExcelFileName = paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_silhouettes.xlsx"), SheetNames = paste(c('Sk', method, 'sil'), collapse = '_'), AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = NA, FreezeRow = 1, FreezeCol = 1)
  crit.values <- vapply(1:length(k.test), function(k) { return(skmeans.res[[k]]$value) }, .1)
  best.idx <- which.max(crit.values)
  bestK <- k.test[best.idx]
  
  ## Ploting criterion values
  svg(filename = paste0(outdir, '/skmeans_k', min(k.test), '.k', max(k.test), '_method.', methodword, '_maxIter.', maxIter, '_best.', bestK, '.svg'), width = 1024/96, height = 1024/96)
  plot(k.test, crit.values, type = "b", pch = 20, xlab = "K", ylab = "Criterion", main = paste0("skmeans criterion (to maximize)\nmethod = ", methodword))
  points(bestK, crit.values[best.idx], pch = 18, cex = 2, col = 2)
  segments(bestK, 0, bestK, crit.values[best.idx], lty = 2, col = 2)
  # dev.off()
  svg_off()
  
  ## Membership
  skmeans.membership <- data.frame(Sample = colnames(data), foreach(k = k.test, .combine = "cbind") %do% skmeans.res[[as.character(k)]]$cluster, stringsAsFactors = FALSE)
  colnames(skmeans.membership) <- c("Sample", k.test)
  memb_root <- paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_membership")
  memb_root_best <- paste0(memb_root, "_best.", bestK)
  memb_root_all <- paste0(memb_root, "_AllK")
  saveRDS(object = skmeans.membership[, c(1, which(colnames(skmeans.membership) == bestK))], file = paste0(memb_root_best, '.RDS'))
  skm_df <- skmeans.membership[,c(1,which(colnames(skmeans.membership) == bestK))]
  write.table(x = skm_df, file = paste0(memb_root_best, '.txt'), sep = "\t", quote = FALSE, row.names = FALSE)
  WriteXLS::WriteXLS(x = skm_df, ExcelFileName = paste0(memb_root_best, '.xlsx'), SheetNames = paste('Membership'), AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = NA, FreezeRow = 1, FreezeCol = 1)
  write.table(x = skmeans.membership, file = paste0(memb_root_all, '.txt'), sep = "\t", quote = FALSE, row.names = FALSE)
  WriteXLS::WriteXLS(x = skmeans.membership, ExcelFileName = paste0(memb_root_all, '.xlsx'), SheetNames = paste('Membership'), AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = NA, FreezeRow = 1, FreezeCol = 1)
  
  ## AVG plot
  bestk.class <- skmeans.membership[,which(colnames(skmeans.membership) == bestK)]
  medsampz <- foreach(k = unique(bestk.class), .combine = "cbind") %do% { return(vapply(1:nrow(data), function(x) { return(median(data[x,bestk.class == k])) }, .1)) }
  # png(paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_best.", bestK, "_mediansamples.png"), width=1650, 1024)
  svg(paste0(memb_root_best, '_mediansamples.svg'), width = 1850/96, height = 1024/96)
  plot(0,0, type = "n", xlim = c(1, nrow(medsampz)), ylim = range(medsampz), xaxs = "i", xlab = "Value index", ylab = "Value", main = paste0("SKmeans (", method, ") median samples for best K (", bestK, ").\nPopulations = ", paste0(as.vector(table(bestk.class)), collapse = ", ")))
  for (k in 1:ncol(medsampz)) lines(medsampz[,k], col = k, lwd = 2)
  # dev.off()
  svg_off()
  
  # msrange <- vapply(1:nrow(medsampz), function(x) {diff(range(medsampz[x,]))}, .1)
  # png(paste0(outdir, "/skmeans_k", min(k.test), ".k", max(k.test), "_method.", methodword, "_maxIter.", maxIter, "_best.", bestK, "_mediansamples_sorted.png"), width=1650, 1024)
  # plot(0,0, type = "n", xlim = c(1,nrow(medsampz)), ylim = range(medsampz), xaxs = "i", xlab = "Value index", ylab = "Value", main = paste0("SKmeans (", method, ") median samples for best K (", bestK, ").\nPopulations = ", paste0(as.vector(table(bestk.class)), collapse = ", ")))
  # plot(0,0, type = "n", xlim = c(1,100), ylim = range(medsampz), xaxs = "i", xlab = "Value index", ylab = "Value", main = paste0("SKmeans (", method, ") median samples for best K (", bestK, ").\nPopulations = ", paste0(as.vector(table(bestk.class)), collapse = ", ")))
  # for (k in 1:ncol(medsampz)) lines(medsampz[order(msrange, decreasing = TRUE),k], col = k)
  # dev.off()
}

## CONTROL examples
### method = 'genetic' (default)
# control <- list(maxiter = 12, popsize = 6, mutations = .1)
### method = 'pclust'
# control <- list(maxiter = 100, nruns = 1, maxchains = 0)
### method = 'CLUTO'
### CLUTO methods : rbr, rb, direct, agglo, graph, bagglo
# control <- list(vcluster = '/path/to/vcluster', colmodel = '', control = paste0('-clmethod=', clutomethod))
### method = 'gmeans'
### gmneans methods :
### s: spherical k-means algorithm (default)
### e: euclidean k-means algorithm
### b: information bottleneck algorithm
### k: kullback_leibler k-means algorithm
### d: diametric k-means algorithm
# control <- list(gmeans = '/path/to/gmeans', control = paste0(' -a ', gmeans.a))
### method = 'kmndirs'
# control <- list(nstart = 10, maxiter = 10)


# source("/home/job/svn/genomics/CGH/R/chrload.R")

# k.test <- 2:15
# maxIter <- 1000
# method <- NULL ### NULL corresponds to standard spherical kmeans. Other available are : genetic, pclust
# method <- "genetic"
# method <- "pclust"
# method <- "kmndirs"
# method <- "gmeans"; gmeans <- "/home/job/Tools/gmeans/gmeans-"; gmeans.a <- "s"
# method <- "CLUTO"; vcluster <- "/home/job/Tools/CLUTO/cluto-2.1.2/Linux-x86_64/vcluster"; clmethod <- "rbr"
# method <- "CLUTO"; vcluster <- "/home/job/Tools/CLUTO/cluto-2.1.2/Linux-x86_64/vcluster"; clmethod <- "rb"
# method <- "CLUTO"; vcluster <- "/home/job/Tools/CLUTO/cluto-2.1.2/Linux-x86_64/vcluster"; clmethod <- "direct"
# method <- "CLUTO"; vcluster <- "/home/job/Tools/CLUTO/cluto-2.1.2/Linux-x86_64/vcluster"; clmethod <- "agglo"
# method <- "CLUTO"; vcluster <- "/home/job/Tools/CLUTO/cluto-2.1.2/Linux-x86_64/vcluster"; clmethod <- "graph"
# method <- "CLUTO"; vcluster <- "/home/job/Tools/CLUTO/cluto-2.1.2/Linux-x86_64/vcluster"; clmethod <- "bagglo"

## RUN

# mydatfile <- "/mnt/data_cigogne/job/cna_joint_admixed_analysis/04_most_centered_peak_centering/multipcf_ADMIX2_sub_GC100000_l2r_125s/OS2K6_ADMIX2_sub_Gpeak18_mcpc_125s_hg19_20160621/Data/OS2K6_ADMIX2_Gpeak_mcpc_125s_hg19_20160621.reg"
# mydf <- read.table(mydatfile, header = TRUE, sep="\t", check.names = FALSE, as.is = TRUE)
# mymat <- as.matrix(mydf[,-c(1:9)])


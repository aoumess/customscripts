## NOTE (20241003)
## =============== ##
## Due to a grid bug in latest stable CRAN version (0.28), installation from github is recommended :
## remotes::install_github("renozao/pkgmaker")
## remotes::install_github("renozao/NMF", ref = "devel")

## Other packages to test : aged, NNLM, RcppML (or via scatter::calculateNMF), dimRed

## NMF_run
##
## DESCRIPTION : A wrapper to ease the use of the NMF package for R. This allows
##               NMF clustering using various methods, along with capturing the
##               features/variables contributing to the different clusters, and
##               assess the clusters to a null distribution by shuffling data.
##
## AUTHOR:
## Bastien JOB (bastien.job@gustaveroussy.fr)
##
## DEPENDS ON:
## . NMF
## . foreach
##
## SUCCESSFULLY TESTED R VERSIONS:
## . 3.2.0 to 3.3.0
##
## VERSION NOTES:
##
## v2.0 20241005
## . Almost complete rewrite with updated code to fit most recent version of NMF
## . Output structure (both data objects and data on disc) reformatted
## . Added svg_png.R support for graphical ouptuts that can be edited.
## . Removed the method evaluation part (not needed)
##
## v1.1b 20160621
## . Added a parameter to modifiy the nmf.option("shared.memory") to inactivate shared
##   memory through big.matrix as it makes nmf crash on R3.3.0 under Ubuntu 14.04 with
##   NMF package version 0.20.6.
## . Also added a parameter to set the seed method.
## . Changed the default seed method used to "none", so that the default fixed seed
##   123456 is actually used (was not the case before this!)
##
## v1.1 20160426
## . Limited the methods for comparison to six methods (others made comparison fail)
## . Changed default maxIter to 5000 (from 3000)
## . Added a help function
## . Corrected a bug where features scores and contributions to metagenes were not saved

## v1.0 20160302
## . First release

source("/home/job/gits/customscripts/R/svg_png.R")

## NMF methods list ====
## Just gives the names of built-in NMF methods in package NMF
## NOTA : "ls-nmf" and "pe-nmf" seem to fail with my typical datasets...
nmf_methods <- function(echo = FALSE) {
  nmfm <- NMF::nmfAlgorithm()
  if (echo) print(nmfm)
  return(nmfm)
}

## Run NMF ====
nmf_run <- function(data = NULL, ranks = 2:3, method = "brunet", default_seed_method = "none", my_seed = 1337, nrun = 100, maxIter = 5000, shift.if.neg = TRUE, classes = NULL, shuffle = TRUE, shuffle.factor = .25, ncores = NULL, odir = getwd(), shared.memory = TRUE) {
  if (!is.matrix(data)) stop("'data' must be a numeric matrix !")
  if (!method %in% NMF::nmfAlgorithm()) stop("Given NMF algorithm is not supported ! See nmf_methods()")
  if (!is.null(classes)) {
    if (nrow(classes) != ncol(data)) stop("'classes' length and number of columns in 'data' do not match !")
    for (x in colnames(classes)) classes[[x]] <- as.factor(classes[[x]])
  }
  if (length(ranks) < 2) stop('At least two rank values shoud be given !')
  if (any(ranks < 2)) stop('Ranks should start at 2 at least')
  if (!all(is.integer(ranks))) stop('All ranks should be integers !')
  if (!is.null(ncores) & !is.numeric(ncores)) stop("'ncores' must be a positive numeric value !")
  if (length(method) != 1) stop('A single NMF method name should be given !')
  if (!dir.exists(odir)) stop("'odir' must be an existing and accessible directory !")
  datmin <- min(data, na.rm = TRUE)
  if (shift.if.neg & datmin < 0) data <- data - datmin
  
  ## Set options
  if (!is.null(ncores)) NMF::nmf.options("cores" = ncores)
  NMF::nmf.options("maxIter" = maxIter)
  NMF::nmf.options("random.seed" = my_seed)
  NMF::nmf.options("grid.patch" = TRUE) ## To avoid PDF blank page
  NMF::nmf.options("shared.memory" = shared.memory)
  
  ## Create output dir
  methodword <- sub(pattern = "/", replacement = "-", x = method)
  oroot <- paste(c("NMF", paste0(nrow(data), ".", ncol(data)), methodword, paste0("r", paste(range(ranks), collapse = ".")), paste0("n", nrun), paste0("mI", maxIter)), collapse = '_')
  odir <- paste0(odir, "/", oroot, "/")
  dir.create(path = odir, recursive = TRUE)
  
  ## NMF
  message("Running NMF ...")
  library(NMF)
  nmf.res <- NMF::nmf(x = data, rank = ranks, method = method, seed = my_seed, nrun = nrun)
  
  ## Add data
  nmf.res$data <- data
  ## Add samples membership
  nmf.res$membership <- data.frame( Samples = colnames(nmf.res$fit[[1]]@fit@H), lapply(seq_along(nmf.res$measures$rank), function(r) { as.numeric(NMF::predict(nmf.res$fit[[r]])) }) , stringsAsFactors = FALSE)
  colnames(nmf.res$membership) <- c('Sample', paste0('rank.', ranks))
  ## Add feature scores
  nmf.res$feature.score <- lapply(seq_along(nmf.res$fit), function(ra) { NMF::featureScore(object = nmf.res$fit[[ra]], method = 'kim') })
  names(nmf.res$feature.score) <- ranks
  ## Add selected features
  nmf.res$feature.selected <- lapply(seq_along(nmf.res$fit), function(ra) { NMF::extractFeatures(object = nmf.res$fit[[ra]], method = 'kim') })
  for (x in seq_along(ranks)) names(nmf.res$feature.selected[[x]]) <- seq.int(1:ranks[x])
  ## Add options
  nmf.res$options <- NMF::nmf.options()
  
  ## Save the NMF object
  message("Saving results ...")
  saveRDS(object = nmf.res, file = paste0(odir, '/', oroot, '.RDS'), compress = 'bzip2')
  
  ## Saving clustering metrics, membership
  write.table(nmf.res$measures, paste0(odir, '/', oroot, '_measures.tsv'), sep="\t", quote = FALSE, row.names = FALSE)
  WriteXLS::WriteXLS(x = nmf.res$measures, ExcelFileName = paste0(odir, '/', oroot, '_measures.xlsx'), SheetNames = 'NMF.measures', AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = NA, FreezeRow = 1, FreezeCol = 1)
  write.table(x = nmf.res$membership, file = paste0(odir, '/', oroot, '_membership.tsv'), sep = "\t", quote = FALSE, row.names = FALSE)
  WriteXLS::WriteXLS(x = nmf.res$membership, ExcelFileName = paste0(odir, '/', oroot, '_membership.xlsx'), SheetNames = 'NMF.membership', AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, na = NA, FreezeRow = 1, FreezeCol = 1)
  
  ## Consensus heatmap
  message("Plotting consensus heatmap ...")
  svg(filename = paste0(odir, '/', oroot, '_consensus.hmap.svg'), width = 40, height = 30)
  if (is.null(classes)) NMF::consensusmap(object = nmf.res) else NMF::consensusmap(object = nmf.res, annCol = classes)
  svg_off()
  
  ## Rank plots
  ### Rank metrics
  message("Plotting ranks measures ...")
  svg(filename = paste0(odir, '/', oroot, '_ranks.survey.svg'), width = 15, height = 10)
  rsplot <- plot(nmf.res, what = c("all", "cophenetic", "rss", "residuals", "dispersion", "evar", "sparseness", "sparseness.basis", "sparseness.coef", "silhouette", "silhouette.coef", "silhouette.basis", "silhouette.consensus"))
  plot(rsplot)
  svg_off()
  
  ## Per rank
  message("Generating rank_level output ...")
  for (ar in seq_along(ranks)) {
    crank <- ranks[ar]
    rank_dir <- paste0(odir, '/rank.', crank)
    dir.create(path = rank_dir, recursive = TRUE)
    ## Consensus heatmap
    svg(filename = paste0(rank_dir, '/rank.', crank, '_consensus.hmap.svg'), width = 22, height = 20)
    if (is.null(classes)) NMF::consensusmap(object = nmf.res$fit[[ar]]) else NMF::consensusmap(object = nmf.res$fit[[ar]], annCol = classes)
    svg_off()
    ## Mixture heatmap
    svg(filename = paste0(rank_dir, '/rank.', crank, '_mixt.coeff.hmap.svg'), width = 15, height = 10)
    if (is.null(classes)) NMF::coefmap(nmf.res$fit[[ar]], main = paste0('Mixture coefficients\n(rank ', crank, ')')) else NMF::coefmap(nmf.res$fit[[ar]], main = paste0('Mixture coefficients\n(rank ', crank, ')'), annCol = classes)
    svg_off()
    ## Basis map
    my_dist <- 'pearson'
    my_clust <- 'ward.D'
    svg(filename = paste0(rank_dir, '/rank.', crank, '_mixt.basis.hmap_', my_dist, '.', my_clust, '_all.svg'), width = 15, height = 22)
    NMF::basismap(nmf.res$fit[[ar]], subsetRow = FALSE, distfun = my_dist, hclustfun = my_clust, treeheight = 0, annRow = list(Metagene=':basis'), main = paste0('Basis components for ALL features\n(rank ', crank, ')'))
    svg_off()
    svg(filename = paste0(rank_dir, '/rank.', crank, '_mixt.basis.hmap_', my_dist, '.', my_clust, '_selected.svg'), width = 10, height = 15)
    NMF::basismap(nmf.res$fit[[ar]], subsetRow = TRUE, distfun = my_dist, hclustfun = my_clust, treeheight = 0, annRow = list(Metagene=':basis'), main = paste0('Basis components for SELECTED features\n(rank ', crank, ')'))
    svg_off()
    ## Save selected features
    bf_df <- data.frame(Feature = unname(unlist(nmf.res$feature.selected[[ar]])), Cluster = rep.int(x = 1:crank, times = sapply(nmf.res$feature.selected[[ar]], length)))
    bf_df <- bf_df[!is.na(bf_df$Feature),]
    bf_df$Feature <- rownames(nmf.res$data)[bf_df$Feature]
    WriteXLS::WriteXLS(x = bf_df, ExcelFileName = paste0(rank_dir, '/rank.', crank, '_selected.features.xlsx'), AdjWidth = TRUE, BoldHeaderRow = TRUE)
    ## Save membership
    memb_df <- nmf.res$membership[,c(1,ar)]
    WriteXLS::WriteXLS(x = memb_df, ExcelFileName = paste0(rank_dir, '/rank.', crank, '_samples.membership.xlsx'), AdjWidth = TRUE, BoldHeaderRow = TRUE)
  }
  
  ## Do shuffling ?
  if (shuffle) {
    ## NMF shuffled
    message("Re-running with shuffled data ...")
    set.seed(my_seed)
    data.rand <- NMF::randomize(data)
    nmf.res.rand <- NMF::nmf(x = data.rand, rank = ranks, method = method, seed = my_seed, nrun = round(nrun * shuffle.factor))
    
    ## Save results
    message("Saving shuffled results ...")
    saveRDS(object = nmf.res.rand, file = paste0(odir, '/', oroot, '_SHUFFLED.RDS'), compress = 'bzip2')
    
    ## Ranks plot
    message("Plotting shuffled results ...")
    svg(paste0(odir, '/', oroot, '_ranks.survey_SHUFFLED.svg'), width = 15, height = 10)
    rscplot <- plot(nmf.res, nmf.res.rand)
    plot(rscplot)
    svg_off()
  }
}


## HELP ====
## Help for main function
help.nmf_run <- function() {
  cat("
nmf_run <- function(data = NULL, ranks = 2:3, method = \"brunet\", default_seed_method = \"none\", seed = 123456, nrun = 100, maxIter = 5000, shift.if.neg = TRUE, classes = NULL, shuffle = TRUE, ncores = NULL, odir = NULL, shared.memory = TRUE)

data            :  (matrix)             A features (rows) by samples (cols) matrix of numerical values.
ranks            :  (integer)            A vector of integers corresponding to desired clusters.
method          :  (character)          The NMF method to use. Run nmf_methods() to obtain a list of available methods.
default_seed_method    :  (character)          The default seed method.
seed            :  (numeric)            The random seed.
nrun            :  (integer)            Nomber of runs to proceed to compute a consensus.
maxIter         :  (integer)            Number of maximum iterations for convergence.
shift.if.neg    :  (logical)            Perform a shift of the matrix if negative values are found.
classes         :  (character/factor)   Known classes for the assessed samples, to compare with identified clusters.
shuffle         :  (logical)            Rerun NMF with shuffled data to compare robustness against a null distribution of your data.
suffle.factor   :  (0<numeric<=1)       The proportion of nrun to use for shuffling (only used when shuffle is TRUE)
ncores          :  (integer)            Number of threads to use (by default, all available minus one).
odir            :  (character)          The directory where results will be put.
shared.memory   :  (logical)            Activate shared memory (using big.matrix). Set it to FALSE in case of crash related to this issue.

")
}




#######
## DEMO FROM THE NMF VIGNETTE (old) ====

# data("esGolub")
# 
# 
# ## Reduce (subseting)
# # esGolub <- esGolub[1:200, ]
# 
# ## Remove the unwanted Sample variable
# esGolub$Sample <- NULL
# 
# ## Single run (default)
# system.time(res <- nmf(esGolub, 3))
# 
# ## Retrieving of the fitted model
# fit(res)
# 
# ## Retrieving of the target matrix
# V.hat <- fitted(res)
# dim(V.hat)
# 
# ## Retrieving of the basis matrix (W, giving the metagenes)
# w <- basis(res)
# dim(w)
# 
# ## Retrieving of the mixture coefficient matrix (H, giving the metagenes)
# h <- coef(res)
# dim(h)
# 
# ## Quality and performance
# summary(res)
# summary(res, target = esGolub)
# 
# ## Getting genes scores (contribution to the metagenes)
# s <- featureScore(res)
# 
# ## Getting genes classified to their metagene
# ## NOTA : Only selecting contributing genes !
# sclass <- extractFeatures(res)
# 
# 
# 
# ## SEEDING METHODS (for random seeding)
# ## NOTA : According to the vignette, if the seeding method is deterministic, no need for multiple NMF runs !!
# nmfSeed()
# 
# 
# ## Multiple ranks testing
# ncores <- 6
# estim.r <- nmf(x = esGolub, rank = 2:6, nrun = 10, seed = 123456, .opt = paste0("vp", ncores))
# plot(estim.r)
# consensusmap(estim.r)
# V.random <- randomize(esGolub)
# system.time(estim.r.rand <- nmf(x = V.random, rank = 2:6, nrun = 10, seed = 123456, .opt = paste0("vp", ncores)))
# plot(estim.r, estim.r.rand)
# 
# # estim.mono <- nmf(x = esGolub, rank = 2:6, nrun = 1, seed = 123456, .opt = paste0("vp", ncores))
#   
# ## Comparison of different algorithms (works for a single rank)
# res.multi.method <- nmf(x = esGolub, rank = 3, nrun = 10, list("brunet", "lee", "ns"), seed = 123456, .opt = paste0("vp", ncores))

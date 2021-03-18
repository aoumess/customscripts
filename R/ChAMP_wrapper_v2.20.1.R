###################
## ChAMP_wrapper ##
###################
## A R wrapper to use ChAMP on illumina methylation microarrays (450K or EPIC designs)
## Tested on R v4.0.4 with ChAMP v2.20.1



## FUNCTIONS
############

### Load datafrom IDAT dir
data.loader <- function(in.dir = NULL, out.dir = NULL) {
  ## Checks
  argg <- c(as.list(environment()))
  if(any(is.null(unlist(argg)))) stop('All parameters to the "data.loader" function should be set.')
  library(ChAMP)
  ## Creating output directory if not existing
  dir.create(path = out.dir, recursive = TRUE, showWarnings = FALSE)
  ## Loading data from idat files
  mych3 <- ChAMP::champ.load(directory = in.dir)
  ## Fixing 'Slide' annotation (to check batch effect)
  mych3$pd$Slide <- as.factor(as.character(mych3$pd$Slide))
  ## Saving loaded data
  saveRDS(mych3, file = paste0(out.dir, "/champ.load.RDS"), compress = 'bzip2')
  return(mych3)
}

## Normalization
data.normalizer <- function(raw.data = NULL, norm.method = 'BMIQ', out.dir = NULL) {
  ## Checks
  argg <- c(as.list(environment()))
  if(any(is.null(unlist(argg)))) stop('All parameters to the "data.normalizer" function should be set.')
  ## Creating outpit dirs
  normplot.dir <- paste0(out.dir, '/Normalization_plots')
  dir.create(path = normplot.dir, recursive = TRUE)
  ## Customizing plot options
  if(norm.method == "BMIQ") plot.bmiq <- TRUE else plot.bmiq <- FALSE
  ## Perform normalization
  mych3.norm <- champ.norm(beta = mych3$beta, rgSet = mych3$rgSet, mset = mych3$mset, method = norm.method, resultsDir = normplot.dir, plotBMIQ = plot.bmiq, cores = nthread)
  ## Saving normalized object
  saveRDS(mych3.norm, file = paste0(out.dir, "/champ.norm.RDS"), compress = 'bzip2')
  return(mych3.norm)
}

### Plot SVD for bias factor
svd.plotter <- function(raw.data = NULL, norm.data = NULL, out.dir = NULL) {
  ## Checks
  argg <- c(as.list(environment()))
  if(any(is.null(unlist(argg)))) stop('All parameters to the "svd.plotter" function should be set.')
  dir.create(path = out.dir, recursive = TRUE, showWarnings = FALSE)
  png(paste0(out.dir, "/SVD.png"), 1024, 512)
  cSVD.res <- champ.SVD                        (beta = norm.data, rgSet = raw.data$rgSet, pd = raw.data$pd, resultsDir = out.dir, RGEffect = TRUE, PDFplot = FALSE, Rplot = TRUE)
  dev.off()
  return(invisible(NULL))
}

### Perform QC plots
qc.plotter <- function(out.dir = NULL, raw.data = NULL, norm.data = NULL, pheno = NULL) {
  ## Checks
  argg <- c(as.list(environment()))
  if(any(is.null(unlist(argg)))) stop('All parameters to the "qc.plotter" function should be set.')
  # if(any(is.null(c(out.dir, raw.data, norm.data, pheno)))) stop('All parameters to the "qc.plotter" function should be set.')
  ## Creating raw qc output directory
  qc.raw.dir <- paste0(out.dir, "/QC_plots/raw")
  dir.create(qc.raw.dir, recursive = TRUE)
  ## Perform QC on raw data
  ChAMP::champ.QC(beta = raw.data, pheno = pheno, resultsDir = qc.raw.dir, Feature.sel = "None", Rplot = FALSE)
  ## Creating normalized qc output directory
  qc.norm.dir <- paste0(out.dir, "/QC_plots/normalized")
  dir.create(qc.norm.dir, recursive = TRUE)
  ## Perform QC on normalized data
  ChAMP::champ.QC(beta = norm.data, pheno = pheno, resultsDir = qc.norm.dir, Feature.sel = "None", Rplot = FALSE)
  return(invisible(NULL))
}

### Perform GSEA analysis (only tested with method = 'ebayes')
### '...' = any additional parameter to ChAMP::champ.GSEA(), namely 'DMP' or 'DMR', depending on the type of input object. See ?ChAMP::champ.GSEA() and its examples.
GSEA.runner <- function(beta = NULL, pheno = NULL, method = 'ebayes', adjPval = .05, out.dir = NULL, file.rootname = 'GSEA', cores = 1, ...) {
  ## Checks
  argg <- c(as.list(environment()), list(...))
  if(any(is.null(unlist(argg)))) stop('All parameters to the "GSEA.runner" function should be set.')
  # if(any(is.null(c(beta, pheno, method, adjPval, out.dir, file.rootname, cores))))
  ## Perform GSEA
  ebGSEA.RES <- ChAMP::champ.GSEA(beta = beta, pheno = pheno, method = method, adjPval = adjPval, cores = cores, ...)
  out.root <- paste0(out.dir, '/', paste(c(file.rootname, 'GSEA', method), collapse = '_'))
  ## Save results object
  saveRDS(ebGSEA.RES, file = paste0(out.root, '.RDS'), compress = 'bzip2')
  ## If there are some significant results
  if (length(ebGSEA.RES$GSEA) > 0) {
    ebGSEA.RES.df <- ebGSEA.RES$GSEA[['Rank(P)']]
    colnames(ebGSEA.RES.df) <- c('Enriched_Genes_Nbr', 'Wilcoxon_AUC', 'Wilcoxon_RawP', 'KPMT_RawP', 'AdjP.BH')
    ## Output significant results as a table
    write.table(x = data.frame(TERM = rownames(ebGSEA.RES$GSEA[['Rank(P)']]), ebGSEA.RES.df, ENR.GENES = vapply(rownames(ebGSEA.RES.df), function(x) { paste(ebGSEA.RES$EnrichGene[[x]], collapse = ',')}, 'a')), file = paste0(out.root, '_TERMS_results.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
    ## Output gene-level results
    write.table(x = data.frame(SYMBOL = rownames(ebGSEA.RES$gtResult), ebGSEA.RES$gtResult), file = paste0(out.root, '_GENES_results.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
  }
  return(invisible(NULL))
}




## VARIABLES
idat.dir <- '/home/job/WORKSPACE/TNE_Ivana/methylation_450k/IDAT_192/'
res.dir <- '/home/job/WORKSPACE/TNE_Ivana/methylation_450k/ChAMP/NORMAL_48/'
norm.method <- 'BMIQ'
gsea.method <- 'ebayes'
nthread <- 2
min.pval <- 5E-02
adjust.method <- 'BH'
diff.categories <- c('Histology', 'Gender')
cna.do <- FALSE
dma.do <- TRUE
qcplot.do <- TRUE

## Getting arguments just to know if at least one was supplied
cmd.args <- base::commandArgs(trailingOnly=TRUE)

## Actually parse arguments
option_list <- list(
  optparse::make_option(c('-i', '--idat-dir'), type = 'character', default = getwd(), help = "Directory containing idat files (and samplesheet) [default = '%default']", metavar = "[PATH]"),
  optparse::make_option(c('-o', '--out-dir'), type = 'character', default = getwd(), help = "Output file name [default = '%default']", metavar = '[PATH]'),
  optparse::make_option(c('-n', '--norm-method'), type = 'character', default = 'BMIQ', help = "Normalization method. Should be one of 'BMIQ', 'SWAN', 'PBC', 'FunctionalNormalize' [default = '%default']", metavar = '[NAME]'),
  optparse::make_option(c('-d', '--dma'), type = 'logical', default = TRUE, help = 'Perform differential methylation analysis (DMA) [default = %default]', metavar = '[TRUE|FALSE]'),
  optparse::make_option(c('-c', '--category-names'), type = 'character', default = 'Slide', help = 'Character, or vector of characters (separated by commas) corresponding to the column names in the samplesheet to perform differential analysis', metavar = '[column_names]'),
  optparse::make_option(c('-g', '--gsea'), type = 'logical', default = TRUE, help = 'Perform gene-set enrichment analysis (GSEA) on DMA results [default = %default]', metavar = '[TRUE|FALSE]'),
  optparse::make_option('--gsea-method', type = 'character', default = 'ebayes', help = "GSEA method. Should be one of 'ebayes', 'fisher', 'gometh' [default = '%default']", metavar = '[NAME]'),
  optparse::make_option(c('-p', '--min-adj-p'), type = 'numeric', default = .05, help = 'Cut-off for FDR-adjusted p-values (for both DMA and GSEA) [default = %default]', metavar = '[0<value<=1]'),
  optparse::make_option(c('-f', '--fdr-method'), type = 'character', default = 'BH', help = 'P-values FDR-adjustment method [default = %default]', metavar = '[BH]'),
  optparse::make_option('--cna', type = 'logical', default = TRUE, help = 'Perform CNA analysis [default = %default]', metavar = '[TRUE|FALSE]'),
  optparse::make_option(c('-q', '--qc-plots'), type = 'logical', default = TRUE, help = 'Perform QC plots [default = %default]', metavar = '[TRUE|FALSE]'),
  optparse::make_option(c('-t', '--thread'), type = 'integer', default = 1, help = 'Number of threads for parallelization. Remember that each thread needs up to 12 GB of RAM for 450k arrays! [default = %default]', metavar = '[INTEGER]')
)

opt_parser <- optparse::OptionParser(option_list = option_list)
suppressWarnings(opt <- optparse::parse_args(opt_parser))

## RUN
######

## If no argument was given, stop with help message
if(length(cmd.args) == 0) {
  optparse::print_help(opt_parser)
} else{
  ## else, let it gooooooo !
  
  ## Fixing category cols
  opt[['category-names']] <- unlist(strsplit(x = as.character(opt[['category-cols']]), ','))
  
  library(ChAMP)
  
  ## Converting opt variables to script variables
  res.dir <- opt[['out-dir']]
  idat.dir <- opt[['idat-dir']]
  norm.method <- opt[['norm-method']]
  gsea.method <- opt[['gsea-method']]
  nthread <- opt[['thread']]
  min.pval <- opt[['min-adj-p']]
  adjust.method <- opt[['fdr-method']]
  diff.categories <-opt[['category-names']]
  cna.do <- opt[['cna']]
  dma.do <- opt[['dma']]
  qcplot.do <- opt[['qc-plots']]
  
  ## LOADING DATA (158s for TNE-POUMON)
  load.dir <- paste0(res.dir, "/01_data_load")
  mych3 <- data.loader(in.dir = idat.dir, out.dir = load.dir)
  
  ## NORMALIZATION (608s for TNE-POUMON)
  norm.dir <- paste0(res.dir, "/02_norm_", norm.method)
  mych3.norm <- data.normalizer(raw.data = mych3, norm.method = 'BMIQ', out.dir = norm.dir)
  
  ## SVD check
  svd.norm.dir <- paste0(norm.dir, "/SVD")
  svd.plotter(raw.data = mych3, norm.data = mych3.norm, out.dir = svd.norm.dir)
  
  ## DIFFERENTIAL ANALYSES
  if(dma.do & length(diff.categories) > 0) {
    message("Performing DMA ...")
    
    for (my.categ in diff.categories) {
      
      ## Checking if category exists
      if (!my.categ %in% colnames(mych3$pd)) {
        stop(paste0("'", my.categ, "' is not available in the samplesheet !"))
      }
      
      ## Go
      dma.dir <- paste0(norm.dir, '/DMA')
      
      print(paste0("Evaluating [", my.categ, "] ..."))
      
      oripheno <- as.factor(mych3$pd[[my.categ]])
      
      ## Checking if there are at least 2 levels
      if (nlevels(oripheno) > 1) {
        
        levels(oripheno) <- gsub(pattern = "-", replacement = "", x = levels(oripheno), fixed = TRUE)
        nonasamp <- !is.na(oripheno)
        
        ## Creating phenotype-based output dir
        pheno.dir <- paste0(dma.dir, '/', my.categ)
        dir.create(pheno.dir, recursive = TRUE)
        
        ### QC plots
        if(qcplot.do) qc.plotter(out.dir = pheno.dir, raw.data = mych3$beta[,nonasamp], norm.data = mych3.norm[,nonasamp], pheno = oripheno[nonasamp])
        
        ## Generating levels combinations
        pheno.comb <- combn(levels(oripheno), 2)
        
        ## Generating "vs other" combs if more than two levels
        if (nlevels(oripheno) > 2) {
          pheno.comb <- cbind(pheno.comb, matrix(c(rbind(levels(oripheno), rep('All_other', nlevels(oripheno)))), nrow = 2))
        }
        
        ## Looping on combs
        for (myc in seq_len(ncol(pheno.comb))) {
          
          mypheno <- oripheno
          
          ## Checking if we are in a "vs other" comb type
          if (pheno.comb[2,myc] == 'All_other') {
            other.level <- levels(mypheno) != pheno.comb[1,myc]
            levels(mypheno)[other.level] <- 'All_other'
            mypheno <- relevel(mypheno, ref = 'All_other')
          } else {
            levels(mypheno)[!levels(mypheno) %in% pheno.comb[,myc]] <- NA
            mypheno <- relevel(mypheno, ref = pheno.comb[2,myc])
          }
          mypheno.inv <- relevel(mypheno, ref = levels(mypheno)[2])
          
          ## Comparison directory
          comp.name <- paste(pheno.comb[,myc], collapse = '_vs_')
          comp.name.inv <- paste(rev(pheno.comb[,myc]), collapse = '_vs_')
          comp.dir <- paste0(pheno.dir, '/', comp.name)
          dir.create(comp.dir)
          
          ## Getting non-NA values
          nonasamp <- !is.na(mypheno)
          
          ### QC plots
          if(qcplot.do) qc.plotter(out.dir = comp.dir, raw.data = mych3$beta[,nonasamp], norm.data = mych3.norm[,nonasamp], pheno = oripheno[nonasamp])
          
          ## DIFFERENTIAL METHYLATION PROBES (DMP)
          message("Computing DMPs (differentialy methylated probes) ...")

          dmp.dir <- paste0(comp.dir, "/DMProbes")
          dir.create(dmp.dir, recursive = TRUE)

          mych3.norm.DMP <- try(champ.DMP(beta = mych3.norm[,nonasamp], pheno = mypheno[nonasamp], adjPVal = min.pval, adjust.method = adjust.method))

          if(!is(mych3.norm.DMP, class2 = 'try-error')) {
            saveRDS(mych3.norm.DMP, file = paste0(dmp.dir, "/", paste(c(my.categ, comp.name, 'DMP', 'Adjp'), collapse = '_'), min.pval, "_", adjust.method, ".RDS"))
            dmp.df <- data.frame(ProbeName = rownames(mych3.norm.DMP[[1]]), mych3.norm.DMP[[1]])
            write.table(x = dmp.df, file = paste0(dmp.dir, "/", paste(c(my.categ, comp.name, 'DMP', 'Adjp'), collapse = '_'), min.pval, "_", adjust.method, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
            dmp.df$logFC <- -dmp.df$logFC
            dmp.df$t <- -dmp.df$t
            write.table(x = dmp.df, file = paste0(dmp.dir, "/", paste(c(my.categ, comp.name.inv, 'DMP', 'Adjp'), collapse = '_'), min.pval, "_", adjust.method, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
            rm(dmp.df)

            ## GSEA on diff results (using unbiased ebGSEA)
            message("Performing GSEA on DMPs ...")
            GSEA.runner(beta = mych3.norm[, nonasamp], DMP = mych3.norm.DMP[[1]], pheno = mypheno[nonasamp], method = 'ebayes', adjPval = min.pval, out.dir = dmp.dir, file.rootname = paste(c(my.categ, comp.name, 'DMP'), collapse = '_'), cores = nthread)
          }

          rm(mych3.norm.DMP)

          ## DIFFERENTIAL METHYLATION REGIONS (DMR)
          message("Computing DMRs (differentialy methylated regions) ...")
          dmr.bh.dir <- paste0(comp.dir, "/DMRegions/Bumphunter")
          dir.create(dmr.bh.dir, recursive = TRUE)

          ### BUMPHUNTER (takes some time due to 500 bootstraps)
          mych3.norm.DMR.bh <- try(champ.DMR(beta = mych3.norm[,nonasamp], pheno = mypheno[nonasamp], arraytype = "450k", method = "Bumphunter", minProbes = 5, adjPvalDmr = min.pval, B = 500, cores = nthread))

          if(!is(mych3.norm.DMR.bh, class2 = 'try-error')) {
            saveRDS(mych3.norm.DMR.bh, file = paste0(dmr.bh.dir, "/", paste(c(my.categ, comp.name, 'DMR', 'Adjp'), collapse = '_'), min.pval, "_", adjust.method, ".RDS"))
            write.table(x = mych3.norm.DMR.bh, file = paste0(dmr.bh.dir, "/", paste(c(my.categ, comp.name, 'DMR', 'Adjp'), collapse = '_'), min.pval, "_", adjust.method, ".txt"), quote = FALSE, sep ="\t", row.names = FALSE)

            ## GSEA on diff result (using unbiased ebGSEA)
            message("Performing GSEA on DMRs ...")
            GSEA.runner(beta = mych3.norm[, nonasamp], DMR = mych3.norm.DMR.bh, pheno = mypheno[nonasamp], method = gsea.method, adjPval = min.pval, out.dir = dmr.bh.dir, file.rootname = paste(c(my.categ, comp.name, 'DMR'), collapse = '_'), cores = nthread)
          }

          rm(mych3.norm.DMR.bh)

          ## DMR ON INVERTED LEVELS
          ### Gives different results, BUT NO GSEA TO PERFORM, as it outputs the exact same results (as there is no signed score)

          message("Computing DMRs on inverted levels ...")
          mych3.norm.DMR.bh.inv <- try(champ.DMR(beta = mych3.norm[,nonasamp], pheno = mypheno.inv[nonasamp], arraytype = "450k", method = "Bumphunter", minProbes = 5, adjPvalDmr = min.pval, B = 500, cores = nthread))

          if(!is(mych3.norm.DMR.bh.inv, class2 = 'try-error')) {
            saveRDS(mych3.norm.DMR.bh.inv, file = paste0(dmr.bh.dir, "/", paste(c(my.categ, comp.name.inv, 'DMR', 'Adjp'), collapse = '_'), min.pval, "_", adjust.method, ".RDS"))
            write.table(x = mych3.norm.DMR.bh.inv, file = paste0(dmr.bh.dir, "/", paste(c(my.categ, comp.name.inv, 'DMR', 'Adjp'), collapse = '_'), min.pval, "_", adjust.method, ".txt"), quote = FALSE, sep ="\t", row.names = FALSE)
          }

          rm(mych3.norm.DMR.bh.inv)
          
          
          ## DIFFERENTIAL METHYLATION BLOCKS (DMB)
          message("Computing DMBs (differentialy methylated blocks) ...")
          dmb.bh.dir <- paste0(comp.dir, "/DMBlocks/Bumphunter")
          dir.create(dmb.bh.dir, recursive = TRUE)
          
          mych3.norm.DMB.bh <- try(ChAMP::champ.Block(beta = mych3.norm[,nonasamp], pheno = mypheno[nonasamp], arraytype = "450k", B = 500, cores = nthread))
  
          if(!is(mych3.norm.DMB.bh, class2 = 'try-error')) {
            saveRDS(mych3.norm.DMB.bh, file = paste0(dmb.bh.dir, "/", paste(c(my.categ, comp.name, 'DMB', 'Adjp'), collapse = '_'), min.pval, "_", adjust.method, ".RDS"))
            write.table(x = data.frame(Block = rownames(mych3.norm.DMB.bh$Block), mych3.norm.DMB.bh$Block), file = paste0(dmb.bh.dir, "/", paste(c(my.categ, comp.name, 'DMB', 'Adjp'), collapse = '_'), min.pval, ".txt"), quote = FALSE, sep ="\t", row.names = FALSE)
  
            ## NO GSEA ON DMB AS RESULTS ARE IDENTICAL TO DMR
          }
          rm(mych3.norm.DMB.bh)
          
          ## DMB ON INVERTED LEVELS
          ### Gives different results
          message("Computing DMBs on inverted levels ...")
          mych3.norm.DMB.bh.inv <- try(ChAMP::champ.Block(beta = mych3.norm[,nonasamp], pheno = mypheno.inv[nonasamp], arraytype = "450k", B = 500, cores = nthread))
          
          if(!is(mych3.norm.DMB.bh.inv, class2 = 'try-error')) {
            saveRDS(mych3.norm.DMB.bh.inv, file = paste0(dmb.bh.dir, "/", paste(c(my.categ, comp.name.inv, 'DMB', 'Adjp'), collapse = '_'), min.pval, "_", adjust.method, ".RDS"))
            write.table(x = data.frame(Block = rownames(mych3.norm.DMB.bh.inv$Block), mych3.norm.DMB.bh.inv$Block), file = paste0(dmb.bh.dir, "/", paste(c(my.categ, comp.name.inv, 'DMB', 'Adjp'), collapse = '_'), min.pval, ".txt"), quote = FALSE, sep ="\t", row.names = FALSE)
          }
          rm(mych3.norm.DMB.bh.inv)
        }
      }
    }
  }
  
  ## CNA
  if(cna.do) {
    ### Using package internal reference profiles
    cna.def.dir <- paste0(norm.dir, '/CNA/champCtls_Ref')
    dir.create(cna.def.dir)
    system.time(CNA.res.def <- ChAMP::champ.CNA(intensity = mych3$intensity, pheno = mych3$pd$Sample_Group, control = TRUE, controlGroup = "champCtls", sampleCNA = TRUE, groupFreqPlots = FALSE, Rplot = TRUE, PDFplot = TRUE, arraytype = '450K', resultsDir = cna.def.dir))
    saveRDS(CNA.res.def, file = paste0(cna.def.dir, '/CNA_champCtls_results.RDS'), compress = 'bzip2')
    ### Using the mean of all samples as control
    cna.mean.dir <- paste0(norm.dir, '/CNA/Mean_Ref')
    dir.create(cna.mean.dir)
    system.time(CNA.res.mean <- ChAMP::champ.CNA(intensity = mych3$intensity, pheno = mych3$pd$Sample_Group, control = FALSE, sampleCNA = TRUE, groupFreqPlots = FALSE, Rplot = TRUE, PDFplot = TRUE, arraytype = '450K', resultsDir = cna.mean.dir))
    saveRDS(CNA.res.mean, file = paste0(cna.mean.dir, '/CNA_mean_results.RDS'), compress = 'bzip2')
  }
}

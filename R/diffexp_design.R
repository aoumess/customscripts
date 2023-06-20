## This script regroups a set of function to preprocess, analyze and plot RNA expression data, and more especially for gene expression differential analysis (using a design table)

## DEPENDENCES ====
### CRAN
cran_dep <- c('amap', 'circlize', 'coop', 'ggplot2', 'randomcoloR', 'stringr')
### BioConductor
bioc_dep <- c('BiocParallel', 'ComplexHeatmap', 'clusterProfiler', 'DESeq2', 'DOSE', 'enrichplot', 'IHW', 'limma', 'matrixStats', 'MeSHDbi', 'meshes', 'msigdbr', 'ReactomePA', 'sparseMatrixStats', 'SummarizedExperiment')

## Install dependences ====
inst_pack <- utils::installed.packages()[,1]
### CRAN
to_inst <- ! cran_dep %in% inst_pack
if(any(to_inst)) for (x in cran_dep[to_inst]) install.packages(x)
### BiocManager
if(!'BiocManager' %in% inst_pack) install.packages('BiocManager')
to_inst <- ! bioc_dep %in% inst_pack
if(any(to_inst)) for (x in bioc_dep[to_inst]) BiocManager::install(x)



## FUNCTIONS ====

## Create a vector of distinct colors
distinct_color_maker <- function(ncolors = 1, my.seed = 12345) {
  set.seed(my.seed)
  return(randomcoloR::distinctColorPalette(k = ncolors, runTsne = if(ncolors > 50) TRUE else FALSE))
}

### Get medoids (for a PCA plot)
## myvec = values (all)
## splitvec = vector of categories (by ex, a factor)
get_medoid <- function(myvec = NULL, splitvec=NULL) {
  splitvec[is.na(splitvec)] <- 'NA'
  vapply(unique(splitvec), function(x) { median(myvec[splitvec == x], na.rm = TRUE) }, .1)
}

## Drawr boxplots with stripchart, except for outliers [LIST VERSION]
## x is a named list, each element of this list will be a new box;
## col can be a single color (then all boxes will have the same color), or a color vector of length(x)
## line.type can be 'medsd' (then median +/- 1 sd of the global population will be drawn) or 'quartiles' (then median, Q1 & Q3 will be drawn)
boxplot2l <- function(x = list(), col = 'lightgray', fill = NULL, main = 'boxplot2', xlab = '', ylab = 'Y', vertical = TRUE, ylim = NULL, cex = 1, line.type = 'medsd') {
  x.meds <- sapply(x, function(i) { median(i, na.rm = TRUE) })
  glob.med <- median(unlist(x), na.rm = TRUE)
  glob.sd <- sd(unlist(x), na.rm = TRUE)
  x.Q1 <- sapply(x, function(i) { quantile(x = i, probs = .25, na.rm = TRUE) })
  x.Q3 <- sapply(x, function(i) { quantile(x = i, probs = .75, na.rm = TRUE) })
  x.IQR <- x.Q3 - x.Q1
  x.len <- sapply(x, length)
  x.noO <- lapply(seq_along(x), function(i) { x[[i]][x[[i]] >= (x.Q1[i] - (1.5 * x.IQR[i])) &  (x[[i]] <= x.Q3[i] + (1.5 * x.IQR[i]))] })
  names(x.noO) <- names(x)
  ## Setting ylim
  if (is.null(ylim)) ylim = range(unlist(x), na.rm = TRUE)
  marb <- round(max(nchar(names(x)))*.75)
  par(mar=c(marb, 4, 2, 2) + 0.1)
  boxplot(x = x, col = fill, border = col, main = paste0(main, ' (', sum(x.len), ')'), ylab = ylab, names = NA, lwd = 2, ylim = ylim, pch = 20, cex = cex)
  stripchart(x = x.noO, method = 'jitter', vertical = vertical, add = TRUE, col = col, pch = 20, cex = cex)
  if(!is.null(line.type)) { if(line.type == 'medsd') abline(h = glob.med + c(0, glob.sd, -glob.sd), col = 2, lty = c(2,3,3), lwd = 2) else if(line.type == 'quartiles') abline(h = c(x.Q1, x.Q3), col = 2, lty = c(3,3), lwd = 2)
  }
  axis(1, at=1:length(x), labels=paste0(names(x), ' (', x.len, ')'), las = 2)
}

### Convert counts to log10(+1)
## x should be a numeric vector or matrix
int2l10 <- function(x = NULL, epsilon = 1) { log10(x + epsilon) }

## Function to assess the weight of annotation covariates in a (sample x annotations) data.frame, on a (feature x sample) data matrix, through correlation (for a continuous covariate) or Kruskal-Wallis statistic (for factors)
## . mat                [f x s num matrix]      A normalized numeric matrix
## . annot.df           [s x a data.frame]      A data.frame with covariates as columns (numeric or factor)
## . factor.names       [vec(char)]             Column names of annot.df corresponding to factor covariates
## . conti.colnames     [vec(char)]             Column names of annot.df corresponding to contiunous covariates
## . red.method         [char]                  Dimension reduction method ['PCA', 'MDS.euc', 'MDS.spear']
## . ndim.max           [int>0]                 Number of dimensions to compute and plot
## . center             [bool]                  Center mat ?
## . scale              [bool]                  Scale mat ?
## . coef.cut           [0<=num<1]              Do not display coefficients inferior to this value on the heatmap
## . color.palette      [vec(col)]              Color vector (length 2) for the heatmap
## . out.file           [char]                  Output PNG file name (and path)
assess_covar <- function(mat = NULL, annot.df = NULL, factor.names = NULL, conti.names = NULL, red.method = 'pca', ndim.max = 10, center = TRUE, scale = TRUE, coef.cut = 0, color.palette = c("white", "orangered3"), out.file = paste0(getwd(), '/Assess_covariates.png')) {
  ## Checks
  ### Mandatory
  if (is.null(mat)) stop('A (f feature by s sample) matrix [mat] is required.')
  if (is.null(annot.df)) stop('An annotation data.frame [annot.df] is required.')
  if (all(is.null(c(factor.names, conti.names)))) stop('At least one of [factor.names] or [conti.colnames] should not be NULL.')
  if (ndim.max <= 0) stop('[ndim.max] should be a non-null positive integer (and <= s samples).')
  if (!dir.exists(dirname(out.file))) stop('Path to [out.file] does not exist.')
  if (!tolower(red.method) %in% c('pca', 'mds.euc', 'mds.spear')) stop('Unknown reduction method')
  ### Compatibility
  if (ndim.max > ncol(mat)) {
    message('WARNING : requested [ndim.max] is higher than samples in [mat]. Reducing it to [mat] samples.')
    ndim.max <- ncol(mat)-1
  }
  # message(ndim.max)
  if (nrow(annot.df) != ncol(mat)) stop('There should be the same number of samples in [mat] (columns) and [annot.df] (rows)')
  if(!is.null(factor.names)) {
    if(!all(factor.names %in% colnames(annot.df))) stop('All [factor.names] should be in colnames of [annot.df].')
  }
  if(!is.null(conti.names)) {
    if(!all(conti.names %in% colnames(annot.df))) stop('All [conti.names] should be in colnames of [annot.df].')
  }
  ## RUN
  ## Center / scale ?
  if (any(c(center, scale))) mat <- base::scale(x = mat, center = center, scale = scale)
  ## Dimension reduction
  if (tolower(red.method) == 'pca') norm.red <- base::svd(x = mat, nv = ndim.max)$v
  if (tolower(red.method) == 'mds.euc') norm.red <- stats::cmdscale(d = dist(x = t(mat), method = 'euclidean'), k = ndim.max)
  if (tolower(red.method) == 'mds.spear') norm.red <- stats::cmdscale(d = as.dist(1-cor(mat, method = 'spearman')), k = ndim.max)
  col.names <- c(factor.names, conti.names)
  col.types <- c(rep('factor', length(factor.names)), rep('continuous', length(conti.names)))
  ## Setting output matrix
  bc.mat <- matrix(NA, ncol = length(col.names), nrow = ndim.max, dimnames = list(paste0(toupper(red.method), seq_len(ndim.max)), col.names))
  ## Filling matrix
  for (cn in seq_along(col.names)) {
    # message(col.names[cn])
    if (col.names[cn] %in% conti.names) {
      cv2cor <- annot.df[[col.names[cn]]]
      nona <- !is.na(cv2cor)
      bc.mat[, cn] <-  abs(cor(x = cv2cor[nona], y = norm.red[nona,], method = 'spearman'))
    } else if (col.names[cn] %in% factor.names & length(unique(annot.df[[col.names[cn]]])) > 1) {
      b2kw <- annot.df[[col.names[cn]]]
      nona <- !is.na(b2kw)
      for (si in seq_len(ndim.max)) {
        bc.mat[si,cn] <- kruskal.test(x = norm.red[nona,si], g = as.factor(b2kw[nona]))$statistic / nrow(norm.red)
      }
    }
  }
  bc.mat[bc.mat < coef.cut] <- 0
  ## Heatmap
  myRamp.col <- circlize::colorRamp2(c(0, 1), color.palette)
  BC.hm <- ComplexHeatmap::Heatmap(matrix = bc.mat,
                                   name = 'Weight',
                                   col = myRamp.col,
                                   na_col = 'grey75',
                                   cluster_rows = FALSE,
                                   cluster_columns = FALSE,
                                   rect_gp = grid::gpar(col = "darkgrey", lwd=0.5),
                                   column_title = 'Batch factors and covariates weight on dataset',
                                   row_title = paste0(toupper(red.method), ' dimensions'),
                                   column_split = col.types,
                                   top_annotation = ComplexHeatmap::HeatmapAnnotation(Type = col.types, col = list(Type = setNames(object = c('lightblue','pink'), nm = c('factor', 'continuous')))))
  png(filename = out.file, width = 800, heigh = 1000)
  ComplexHeatmap::draw(BC.hm)
  dev.off()
}

## Converts a factor to a design dataframe to be used in DE.test
## Function to generated a "full" comparison design from a minimalistic dataframe containing samplenames and factor(s) to compare, to use in DE.test
## . init_df            [data.frame]    A data.frame with at least 2 columns, one of which should be named as defined in [samples_colname], the others being the factors to use to generate the comparison design
## . samples_colname    [char]          Name of the column to use from init_df as sample names
## . conditions     [vec(char)]         Name of the column(s) to use in the design as conditions to compare to regress. They have to be in [init_df] !
## . covariates     [list(vec(char))]   List of name(s) of the column(s) to use in the design as covariates to regress. They should not be in [init_df] !
## . add_inverted       [logical]         Add the inverted comparisons of the default factor ones (ie, B_vs_A if A_vs_B is the default one)
## . add_others         [logical]         Add the 'X vs Other' type of comparisons if the factor has more than 2 levels
## . only_others        [logical]         Only generate the 'X vs Other' type of comparisons if the factor has more than 2 levels
full_design_generator <- function(init_df = NULL, samples_colname = NULL, conditions = NULL, covariates = NULL, add_inverted = TRUE, add_others = TRUE, only_others = FALSE) {
  
  ## Checks
  if(is.null(init_df)) stop('A starting dataframe containing sample names and factors to compare is  required !')
  if(is.null(samples_colname)) stop('A column name to identify samples in init_df is required !')
  if(!samples_colname %in% colnames(init_df)) stop('Column name [', samples_colname, '] not found in init_df !')
  if(is.null(conditions)) stop('Conditions to compare are required !')
  if (length(conditions) != length(covariates)) stop('Conditions vector and covariates list should have the same length !')
  if(!add_others & only_others) stop("Can't restrict to 'vs_Other' comparisons if add_others is not set to TRUE !")
  if(any(! conditions %in% colnames(init_df))) stop('At least one condition was not found in the header of init_df !')
  # if(any(! unique(unlist(covariates)) %in% colnames(init_df))) stop('At least one covariate was not found in the header of init_df !')
  
  ## Cleaning bad chars
  ### NAMES : columns in initial df
  colnames(init_df) <- gsub(pattern = "\\W", replacement = '.', x = colnames(init_df))
  ### NAMES : samples colname
  samples_colname <- gsub(pattern = "\\W", replacement = '.', x = samples_colname)
  ### NAMES : covariate names
  covariates <- lapply(covariates, function(cova) gsub(pattern = "\\W", replacement = '.', x = cova))
  ### CONTENT : initial df
  for (myc in seq_len(ncol(init_df))) init_df[[myc]] <- gsub(pattern = "\\W", replacement = '.', x = init_df[[myc]])
  
  ### Getting conditions from init df
  factor_colnames <- colnames(init_df)[!colnames(init_df) %in% samples_colname]
  
  des.df <- sapply(seq_along(factor_colnames), function(f) {
    factor_name <- factor_colnames[f]
    message(factor_name)
    my_factor <- init_df[[factor_name]]
    if(!is.factor(my_factor)) my_factor <- as.factor(my_factor)
    
    mylevels <- levels(my_factor)
    combn.res <- combn(rev(mylevels), 2)
    
    if(add_inverted) combn.res <- cbind(combn.res, combn(mylevels, 2))
    all.combz <- sapply(1:ncol(combn.res), function(x) {list(combn.res[1,x], combn.res[2,x])}, simplify = FALSE)
    names(all.combz) <- vapply(1:ncol(combn.res), function(x) { paste(combn.res[, x, drop = TRUE], collapse = '_vs_') }, 'a')
    if(length(mylevels) > 2 & add_others) {
      XvsO.combz <- sapply(as.character(mylevels), function(x) { list(x, mylevels[!mylevels == x])}, simplify = FALSE)
      names(XvsO.combz) <- vapply(mylevels, function(x) { paste0(x, '_vs_Other') }, 'a')
      all.combz <- if(only_others) XvsO.combz else c(all.combz, XvsO.combz)
      rm(XvsO.combz)
    }
  
    ## Filtering comparisons with a class comprised by an unique sample
    for (mycomb in names(all.combz)) {
      class.len <- vapply(all.combz[[mycomb]], function(x) { length(which(my_factor %in% x))}, 1L)
      if(any(class.len == 1)) {
        all.combz[[mycomb]] <- NULL
        message('Removed ', mycomb, ' due to subpopulation(s) of size 1 ...')
      }
    }
    
    ## Generate design table
    f.df <- sapply(seq_along(all.combz), function(x) {
      comb_df <- data.frame(
        'Samples_colname' = samples_colname
        , 'Covar_colnames' = paste(covariates[[f]], collapse = ',')
        , 'Condition_colname' = factor_name
        , 'Condition_A' = paste(all.combz[[x]][[1]], collapse = ',')
        , 'Condition_B' = paste(all.combz[[x]][[2]], collapse = ',')
        , 'Comparison_name' = names(all.combz)[x]
      )
      return(comb_df)
    }, simplify = FALSE)
    f.df = Reduce(f = rbind, x = f.df)
    return(f.df)
  }, simplify = FALSE)
  des.df <- Reduce(f = rbind, x = des.df)
  return(des.df)
}


## FUNCTION TO CONVERT AN msigdbr SPECIES NAME TO AN AnnotationForge ORGANISM PACKAGE ROOTNAME (ie, the package name without the '.db' suffix)
### . 'species' : character ; a species name, as in the 'species_name' column of msigdbr::msigdbr_species()
### NOTE1 : This converter should work with all available species in msigdbr, with the exception of saccharomyces cerevisiae (which does not have 'eg' in its organism name)
### NOTE 2 : See MSigDB collections description here : http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
msigdbr2org <- function(species = NULL) {
  return(paste0('org.', paste0(stringr::str_sub(unlist(base::strsplit(species, ' ')), 1, 1), collapse = ''), '.eg'))
}


## Function to convert any of ENSEMBL / HGNC symbol / ENTREZ to any other
## Returns a df with 3 columns (ENSEMBL, SYMBOL, ENTREZID), with rownames corresponding to the provided features having a correspondance.
ID_converter <- function(features = NULL, feature.type = 'ENSEMBL', species = 'Homo sapiens') {
  exp.types <- c('ENSEMBL', 'SYMBOL', 'ENTREZID')
  feature.type <- toupper(feature.type)
  if (!feature.type %in% exp.types) stop('Feature type should be one of ', paste(exp.types, collapse = ','))
  
  ## Cleaning
  message('Removing duplicated input IDs ...')
  features <- features[!duplicated(features)]
  ## Building the gene ID conversion vector
  toget.types <- exp.types[!exp.types == feature.type]
  gconv <- clusterProfiler::bitr(features, fromType = feature.type, toType = toget.types, OrgDb = paste0(msigdbr2org(species), '.db'))
  ## Strip duplicates
  gconv <- gconv[!duplicated(gconv[[feature.type]]),]
  ## Add features as rownames
  rownames(gconv) <- gconv[[feature.type]]
  ## Add T/F vector
  gconv$Pure.ID <- TRUE
  ## Recall ids from the features vector that are absent from gconv
  misfeat <- features[!features %in% gconv[[feature.type]]]
  misdf <- data.frame(misfeat, misfeat, misfeat, misfeat, FALSE, row.names = 1)
  colnames(misdf) <- colnames(gconv)
  gconv <- rbind(gconv, misdf)
  return(gconv)
}


## Prepare count matrix for DEseq2 test
## exp.mat                matrix(integer)     Sample x feature (gene) raw count matrix. Feature names as rownames.
## dupe.remove            logical             Remove duplicated features (through their rownames)
## novar.filter           logical             Remove features with no variance
## nonzero.minsamp        [0<int<+inf]        Minimum number of samples with non-zero value to keep a feature
DE_prepmatrix <- function(exp.mat = NULL, dupe.remove = TRUE, novar.filter = TRUE, nonzero.minsamp = 5) {
  exp.mat <- round(exp.mat)
  message('Initial sparsity level : ', coop::sparsity(exp.mat))
  message('Initial dimensions : ')
  message(paste(dim(exp.mat), collapse = ', '))
  ## Filter out duped features
  if(dupe.remove) {
    message('Filtering duped features.')
    exp.mat <- exp.mat[!duplicated(rownames(exp.mat)),, drop = FALSE]
    message(paste(dim(exp.mat), collapse = ', '))
    message('\tFiltered sparsity level : ', coop::sparsity(exp.mat))
  }
  ## Filter out features without variance
  if(novar.filter) {
    message('Filtering features without any variance.')
    expsd <- matrixStats::rowSds(exp.mat)
    exp.mat <- exp.mat[!expsd == 0,, drop = FALSE]
    message(paste(dim(exp.mat), collapse = ', '))
    message('\tFiltered sparsity level : ', coop::sparsity(exp.mat))
  }
  ## Filtering features with too many zeros
  message('Filtering features with too many zero count values.')
  nzcount <- ncol(exp.mat) - matrixStats::rowCounts(x = exp.mat, value = 0)
  nzselec <- nzcount >= nonzero.minsamp
  exp.mat <- exp.mat[nzselec,, drop = FALSE]
  message(paste(dim(exp.mat), collapse = ', '))
  message('\tFiltered sparsity level : ', coop::sparsity(exp.mat))
  return(exp.mat)
}

## Perform DE analysis & functional enrichment for contrasts in pairs
## exp.mat                matrix(integer)     Sample x feature (gene) raw count matrix. Feature names as rownames.
## feature.type           character           Type of features given in exp.mat (rownames). Should be one of ['SYMBOL', 'ENSEMBL', 'ENTREZID']. Used to connect features to GSEA/ORA, by converting them to ENTREZID if needed using clusterProfiler::bitr
## annot.df               data.frame          Sample annotations. Should contain a column with the same entries as colnames(exp.mat)
## design.df              data.frame          Design of comparisons to perform. Should contain these column names : [Samples_colname] = Column name of sample names ; [Covar_colnames] = Column name(s) of covariates to regress, coma-separated (can be empty if none to regress); [Condition_colname] = Column name of factor condition to explore for the differential analysis ; [Condition_A] = levels to consider as the condition A (test) ; [Condition_B] = levels to consider as the condition B (ref) ; [Comparison_name] = Name to use for results output.
## assess.factor          logical             Perform assessment of factor covariates using the provided column name(s) corresponding to factor data columns in annot.df
## assess.conti           logical             Perform assessment of continuous covariates using the provided column name(s) corresponding to continuous data columns in annot.df
## adjp.max               0<numeric<1         BH FDR-adjusted p-value cut-off to consider differential genes as significant
## lfc.min                numeric+            Minimal logFoldChange to consider differential genes
## ihw                    logical             Apply Independent Hypothesis Weighting (see IHW::ihw) [TRUE]
## outdir                 character           Path of the output directory
## samples.dist.method    character           Name of the distance method to use for the hierarchical clustering of samples
## samples.hclust.method  character           Name of the aggregation method to use for the hierarchical clustering of samples
## genes.dist.method      character           Name of the distance method to use for the hierarchical clustering of genes
## genes.hclust.method    character           Name of the aggregation method to use for the hierarchical clustering of genes
## msigdb.do              c(bool, bool)       Use the MSigDb dataset collection to perform GSEA/ORA
## go.do                  c(bool, bool)       Use the GeneOntology dataset to perform GSEA/ORA
## do.do                  c(bool, bool)       Use the DiseaseOntology + CancerGeneNetwork + DisGeNet datasets to perform GSEA/ORA
## kegg.do                c(bool, bool)       Use the KEGG dataset to perform GSEA/ORA
## wp.do                  c(bool, bool)       Use the WikiPathways dataset to perform GSEA/ORA
## reactome.do            c(bool, bool)       Use the Reactome dataset to perform GSEA/ORA
## cm.do                  c(bool, bool)       Use the CellMarker dataset to perform GSEA/ORA
## mesh.do                c(bool, bool)       Use the MeSHDb dataset collection to perform GSEA/ORA (warning, this may consume an astronomical amount of RAM. gsea.do is forced to false)
## non.redundant          bool                If any of go.do|kegg.do|wp.do|reactome.do AND msigdb.do are TRUE, discard outputing corresponding versions of the former in MSigDb (as they are corresponding to older, smaller bases)
## species                character           Name of the species analyzed (namely, 'human' or 'mouse')
## enrp.max               0<numeric<1         BH FDR-adjusted p-value cut-off to consider enriched terms as significant
## enr.min.genes          numeric+            Minimum number of significant genes to perform GSEA/ORA
## or.top.max             numeric+            Maximum number of significant genes to consider as a signature for ORA. Also used for plots (boxplots, limited heatmap) using only a portion of all significant genes
## dotplot.maxterms       numeric+            Maximum number of enriched terms to plot in a dotplot (for readability)
## my.seed                numeric             Seed value for RNG (used for GSEA and heatmap annotation colors)
## boxplots               bool                If TRUE, draw boxplots of or.top.max genes
## save.wald              bool                If TRUE, save the DESeq2 object containing the results of the Wald test. This is FALSE by default, as the resulting object can be pretty big.
## color.palette          vec(color)          Vector of 3 colors used for the expression heatmap (lower values, middle, higher)
DE_test <- function(exp.mat = NULL, feature.type = 'SYMBOL', annot.df = NULL, design.df = NULL, assess.factor = NULL, assess.conti = NULL, adjp.max = 5E-02, lfc.min = .7, ihw = TRUE, lfcShrink = TRUE, enrp.max = 1E-02, enr.min.genes = 10, or.top.max = 100, outdir = getwd(), samples.dist.method = 'spearman', samples.hclust.method = 'ward.D', genes.dist.method = 'spearman', genes.hclust.method = 'ward.D', msigdb.do = c(TRUE, TRUE), do.do = c(TRUE, TRUE), go.do = c(TRUE, TRUE), kegg.do = c(TRUE, TRUE), wp.do = c(TRUE, TRUE), reactome.do = c(TRUE, TRUE), mesh.do = c(FALSE, FALSE), non.redundant = TRUE, species = 'Homo sapiens', dotplot.maxterms = 50, my.seed = 1234L, boxplots = TRUE, save.wald = FALSE, heatmap.palette = c("royalblue3", "ivory", "orangered3"), BPPARAM = BiocParallel::SerialParam()) {
  
  if (tolower(species) == 'homo sapiens') {
    Org <- 'org.Hs'
  } else if (tolower(species) == 'mus musculus') {
    Org <- 'org.Mm'
  } else stop("Only 'Homo sapiens' and 'Mus musculus' species are supported !")
  
  ## Cleaning design
  for (x in seq_len(ncol(design.df))) design.df[,x] <- as.character(design.df[,x])
  design.df$Covar_colnames[design.df$Covar_colnames == ''] <- NA
  
  ## CHECKS ====
  ## Mandatory entries
  ### expression matrix
  if (is.null(exp.mat)) stop('An expression matrix is required !')
  if (!is.matrix(exp.mat)) stop('Expression data should be a matrix !')
  ### annotation
  if (is.null(annot.df)) stop('An annotation data.frame is required !')
  if (!is.data.frame(annot.df)) stop('Annotation data should be a data.frame !')
  ### design
  if (is.null(design.df)) stop('A design data.frame is required !')
  if (!is.data.frame(design.df)) stop('Design data should be a data.frame !')
  ## if exp.mat and annot.df have different size
  if (nrow(annot.df) != ncol(exp.mat)) stop("'exp.mat' and 'annot.df' do not have the same amount of samples!")
  ## if some annotation column names provided in design do not exist :
  ### samples :
  if (!all(unique(design.df$Samples_colname) %in% colnames(annot.df))) stop('All provided [Samples_colname] values should be in colnames(annot.df) !')
  ### covariates (outside design):
  if (!all(assess.factor %in% colnames(annot.df))) stop('All provided [assess.factor] values should be in colnames(annot.df) !')
  if (!all(assess.conti %in% colnames(annot.df))) stop('All provided [assess.conti] values should be in colnames(annot.df) !')
  ### covariates (from design):
  covars <- unique(unlist(lapply(seq_len(nrow(design.df)), function(x) unlist(strsplit(x = as.character(design.df$Covar_colnames[x]), split = ',')))))
  covars <- covars[!is.na(covars)]
  if(!all(is.na(covars))) {
    if (!all(covars[!is.na(covars)] %in% colnames(annot.df))) stop('All provided [Covar_colname] values should be in colnames(annot.df) !')
  }
  ### condition :
  if (!all(unique(design.df$Condition_colname) %in% colnames(annot.df))) stop('All provided [Condition_colname] values should be in colnames(annot.df) !')
  ## if samples are not identical
  # exp.mat <- exp.mat[,order(colnames(exp.mat))]
  # annot.df <- annot.df[order(annot.df[[samples.colname]]),]
  # if (!all(colnames(exp.mat) == annot.df[[samples.colname]])) stop(paste0("Content of the sample column '", samples.colname, "' is not identical to 'exp.mat' colnames !"))
  
  ## Get id convertions ====
  sconv <- ID_converter(features = rownames(exp.mat), feature.type = feature.type, species = species)
  exp.types <- c('ENSEMBL', 'SYMBOL', 'ENTREZID')
  other.types <- exp.types[!exp.types == feature.type]
  
  ## Convert main ID to Symbol
  # rownames(exp.mat) <- sconv[rownames(exp.mat), 'SYMBOL']
  
  ## Loading packages ====
  library(DESeq2)
  library(SummarizedExperiment)
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
  
  ## Looping on design entries ====
  for (cur.idx in seq_len(nrow(design.df))) {
    
    ## Limiting annot.df to required columns
    samples.colname <- design.df$Samples_colname[cur.idx]
    cur.cond <- design.df$Condition_colname[cur.idx]
    cur.covars <- unlist(strsplit(x = design.df$Covar_colnames[cur.idx], split = ','))
    cur.covars <- cur.covars[!is.na(cur.covars)]
    cur.condA <- unlist(strsplit(x = design.df$Condition_A[cur.idx], split = ','))
    cur.condB <- unlist(strsplit(x = design.df$Condition_B[cur.idx], split = ','))
    cur.annot.df <- annot.df[, colnames(annot.df) %in% c(samples.colname, cur.covars, cur.cond)]
    cur.name <- design.df$Comparison_name[cur.idx]
    cur.exp.mat <- exp.mat
    cur.annot.df <- annot.df
    
    message('Performing comparison : [', cur.name, '] ...')
    
    ## Filtering special characters in cur.cond
    colnames(cur.annot.df) <- gsub(pattern = "\\W", replacement = '.', x = colnames(cur.annot.df))
    # message(paste(levels(cur.annot.df[[cur.cond]]), collapse = ' ; '))
    levels(cur.annot.df[[cur.cond]]) <- gsub(pattern = "\\W", replacement = '.', x = levels(cur.annot.df[[cur.cond]]))
    # message(paste(levels(cur.annot.df[[cur.cond]]), collapse = ' ; '))
    
    
    ## Adjusting the datasets if needed
    ## Handling NAs in current annotation
    na.check <- is.na(cur.annot.df[[cur.cond]])
    if (any(na.check)) {
      cur.exp.mat <- cur.exp.mat[,!na.check]
      cur.annot.df <- cur.annot.df[!na.check,]
    }
    ## Handling requested levels
    lev.checks <- cur.annot.df[[cur.cond]] %in% c(cur.condA, cur.condB)
    if (!all(lev.checks)) {
      cur.exp.mat <- cur.exp.mat[,lev.checks]
      cur.annot.df <- cur.annot.df[lev.checks,]
    }
    
    ## Forcing a relevel (if some levels were lost)
    # cur.annot.df[[cur.cond]] <- as.factor(as.character(cur.annot.df[[cur.cond]]))
    # levtab <- as.data.frame(table(cur.annot.df[[cur.cond]]), stringsAsFactors = FALSE)
    # levtab <- levtab[order(levtab$Var1),]
    # levtab <- levtab[levtab$Freq > 0,]
    # cur.annot.df[[cur.cond]] <- as.factor(as.character(cur.annot.df[[cur.cond]]))
    # myref <- levels(cur.annot.df[[cur.cond]])[1]
    # cur.annot.df[[cur.cond]] <- relevel(cur.annot.df[[cur.cond]], ref = myref)
    if (!is.factor(cur.annot.df[[cur.cond]])) cur.annot.df[[cur.cond]] <- as.factor(cur.annot.df[[cur.cond]])
    cur.annot.df[[cur.cond]] <- droplevels(cur.annot.df[[cur.cond]])
    
    ## Preparing the table list of all combinations
    # mylevels <- levels(cur.annot.df[[cur.cond]])
    # combn.res <- if(invert.levels) combn(mylevels, 2) else combn(rev(mylevels), 2)
    # all.combz <- sapply(1:ncol(combn.res), function(x) {list(combn.res[1,x], combn.res[2,x])}, simplify = FALSE)
    # names(all.combz) <- vapply(1:ncol(combn.res), function(x) { paste(combn.res[, x, drop = TRUE], collapse = '_vs_') }, 'a')
    # if(length(mylevels) > 2) {
    #   XvsO.combz <- sapply(as.character(mylevels), function(x) { list(x, mylevels[!mylevels == x])}, simplify = FALSE)
    #   names(XvsO.combz) <- vapply(mylevels, function(x) { paste0(x, '_vs_Other') }, 'a')
    #   all.combz <- if(only.others) XvsO.combz else c(all.combz, XvsO.combz)
    #   rm(XvsO.combz)
    # }
    # 
    # ## Filtering comparisons with a class comprised by an unique sample
    # for (mycomb in names(all.combz)) {
    #   class.len <- vapply(all.combz[[mycomb]], function(x) { length(which(cur.annot.df[[cur.cond]] %in% x))}, 1L)
    #   if(any(class.len == 1)) all.combz[[mycomb]] <- NULL
    # }
    
    ## Creating design (factor then Batch, w/o intercept)
    my.textform <- paste(c('~0', unique(c(cur.cond, cur.covars))), collapse = '+')
    
    # my.textform <- paste(c('~0', unique(c(cur.covars, cur.cond))), collapse = '+')
    # my.textform <- paste0('~', paste(unique(c(cur.covars, cur.cond)), collapse = '+'))
    # my.textform <- paste0('~', paste(unique(c(cur.cond, cur.covars)), collapse = '+'))
    # my.textform <- paste(c('~0', unique(c(cur.cond, cur.covars))), collapse = '+')
    my.design <- as.formula(my.textform)
    
    ## Creating the main output dir
    fac.dir <- paste0(outdir, '/Differential_analysis')
    factor.dir <- paste(c(fac.dir, paste0('adjp.', adjp.max, '_lfc.', lfc.min), my.textform), collapse = '/')
    # dir.create(path = factor.dir, recursive = TRUE)
    de.dir <- paste(c(factor.dir, cur.name), collapse = '/')
    dir.create(path = de.dir, recursive = TRUE)
    
    ## Creating the DESeq2 object
    dedf.try <- try(DE2obj <- DESeq2::DESeqDataSetFromMatrix(countData = cur.exp.mat, colData = cur.annot.df, design = my.design))
    if (is(dedf.try, class2 = 'try-error')) {
      message('DEA failed (see DESeq2 error  above) !')
      file.rename(from = de.dir, to = sub(pattern = '/$', replacement = '_FAILED', x = de.dir))
      next()
    }
    rm(cur.exp.mat, cur.annot.df)
    
    ## Saving the DESeq object
    # saveRDS(object = DE2obj, file = paste0(factor.dir, '/', cur.cond, '_rawcounts.RDS'), compress = 'bzip2')
    saveRDS(object = DE2obj, file = paste0(de.dir, '/', cur.cond, '_rawcounts.RDS'), compress = 'bzip2')
    
    ## Normalizing by vst (for PCA & heatmap)
    DE2obj.norm <- DESeq2::vst(object = DE2obj, blind = TRUE)
    norm.mat <- SummarizedExperiment::assay(DE2obj.norm)
    rm(DE2obj.norm)
    
    ## Rough log10 normalization
    # norm.mat <- log10(SummarizedExperiment::assay(DE2obj)+1)
    
    ## Assessing DESIGN covariates, and regressing if requested
    if (length(cur.covars) > 0) {
      ### Assessing covariates
      #### Splitting factor and continuous covariates
      factor.colnames <- conti.colnames <- NULL
      for (cn in cur.covars) if (is.factor(DE2obj@colData[[cn]])) factor.colnames <- c(factor.colnames, cn) else if (is.numeric(DE2obj@colData[[cn]])) conti.colnames <- c(conti.colnames, cn) else stop(paste0('Covariate [', cn, '] is neither a factor nor a numeric/integer vector !'))
      #### Assessing covariates
      try(assess_covar(mat = norm.mat, annot.df = as.data.frame(DE2obj@colData), factor.names = c(cur.cond, factor.colnames), conti.names = conti.colnames, red.method = 'pca', ndim.max = round(ncol(norm.mat)/2), center = TRUE, scale = TRUE, out.file = paste0(de.dir, '/', cur.cond, '_assess_covariates_01_unregressed.png')))
      #### Running limma::removeBatchEffect the good way
      limma.bc.batch2 <- limma.bc.batch1 <- limma.bc.covar <- NULL
      ##### Handling factor covariates
      for (fc in factor.colnames) {
        if (is.null(limma.bc.batch1)) limma.bc.batch1 <- DE2obj@colData[[fc]] else if (is.null(limma.bc.batch2)) limma.bc.batch2 <- DE2obj@colData[[fc]] else message(paste0('Factor [', fc, '] will not be considered for matrix regression by limma::removeBatchEffect as only 2 factors can be used at max.'))
      }
      ##### Handling continuous covariates
      for (cc in conti.colnames) {
        if (is.null(limma.bc.covar)) limma.bc.covar <- as.matrix(DE2obj@colData[, cc, drop = FALSE]) else limma.bc.covar <- cbind(limma.bc.covar, as.matrix(DE2obj@colData[, cc, drop = FALSE]))
      }
      
      ### Performing batch effect regression with limma
      bc.norm.mat <- limma::removeBatchEffect(x = norm.mat, batch = limma.bc.batch1, batch2 = limma.bc.batch2, covariates = limma.bc.covar, design = model.matrix(as.formula(paste0('~0+', cur.cond)), data = DE2obj@colData))
      
      ### Assessing covariates (after regression)
      try(assess_covar(mat = bc.norm.mat, annot.df = as.data.frame(DE2obj@colData), factor.names = c(cur.cond, factor.colnames), conti.names = conti.colnames, red.method = 'pca', ndim.max = round(ncol(norm.mat)/2), center = TRUE, scale = TRUE, out.file = paste0(de.dir, '/', cur.cond, '_assess_covariates_02_regressed.png')))
    }
    
    ### PCAs for DESIGN covariates
    for (p in c(cur.cond, cur.covars)) {
      # png(filename = paste0(factor.dir, '/PCA_vst_', p, '.png'), width = 1100, height = 1000)
      png(filename = paste0(de.dir, '/PCA_vst_', p, '.png'), width = 1100, height = 1000)
      library(ggfortify)
      print(ggplot2::autoplot(prcomp(t(norm.mat)), data = as.data.frame(SummarizedExperiment::colData(DE2obj)), colour = p, size = 3))
      dev.off()
    }
    
    ## Assessing TEST covariates, and regressing if requested
    if (any(!is.null(c(assess.factor, assess.conti)))) {
      # test.mat <- SummarizedExperiment::assay(DE2obj.norm)
      test.mat <- norm.mat
      test.dir <- paste0(de.dir, '/Test_covariates')
      dir.create(test.dir)
      
      tmp.annot <- as.data.frame(SummarizedExperiment::colData(DE2obj))
      ### plot PCA of normalized,unregressed data colored by TEST covariates
      for (p in unique(c(cur.cond, assess.factor, assess.conti))) {
        png(filename = paste0(test.dir, '/PCA_vst_UNREGRESSED_col.', p, '.png'), width = 1100, height = 1000)
        library(ggfortify)
        try(print(ggplot2::autoplot(prcomp(t(test.mat)), data = tmp.annot, colour = p, size = 3)), silent = TRUE)
        dev.off()
      }
      
      ### Plotting UNREGRESSED assessment heatmap
      try(assess_covar(mat = test.mat, annot.df = tmp.annot, factor.names = c(cur.cond, assess.factor), conti.names = assess.conti, red.method = 'pca', ndim.max = round(ncol(test.mat)/2), center = TRUE, scale = TRUE, out.file = paste0(test.dir, '/', cur.cond, '_TEST_covariates_01_unregressed.png')))
      
      #### Handling continuous covariates
      for (cc in assess.conti) {
        message(cc)
        ## Regression the continuous covariate
        ber.try <- try(tmp.mat <- limma::removeBatchEffect(x = test.mat, batch = NULL, batch2 = NULL, covariates = tmp.annot[[cc]], design = model.matrix(as.formula(paste0('~0+', cur.cond)), data = DE2obj@colData)), silent = TRUE)
        if (!is(ber.try, class2 = 'try-error')) {
          ## Plotting REGRESSED assessment heatmap
          try(assess_covar(mat = tmp.mat, annot.df = as.data.frame(DE2obj@colData), factor.names = c(cur.cond, assess.factor), conti.names = assess.conti, red.method = 'pca', ndim.max = round(ncol(tmp.mat)/2), center = TRUE, scale = TRUE, out.file = paste0(test.dir, '/', cur.cond, '_TEST_covariates_02_REGRESSED_', cc, '.png')))
          ### plot PCA of REGRESSED data colored by cur.cond
          png(filename = paste0(test.dir, '/PCA_vst_REGRESSED.', cc, '_col.', cur.cond, '.png'), width = 1100, height = 1000)
          library(ggfortify)
          print(ggplot2::autoplot(prcomp(t(tmp.mat)), data = as.data.frame(SummarizedExperiment::colData(DE2obj)), colour = cur.cond, size = 3))
          dev.off()
        }
      }
      
      #### Handling factor covariates
      for (fc in assess.factor) {
        message(fc)
        ## Regression the continuous covariate
        ber.try <- try(tmp.mat <- limma::removeBatchEffect(x = test.mat, batch = tmp.annot[[fc]], batch2 = NULL, covariates = NULL, design = model.matrix(as.formula(paste0('~0+', cur.cond)), data = DE2obj@colData)), silent = TRUE)
        if (!is(ber.try, class2 = 'try-error')) {
          ## Plotting REGRESSED assessment heatmap
          try(assess_covar(mat = tmp.mat, annot.df = tmp.annot, factor.names = c(cur.cond, assess.factor), conti.names = assess.conti, red.method = 'pca', ndim.max = round(ncol(tmp.mat)/2), center = TRUE, scale = TRUE, out.file = paste0(test.dir, '/', cur.cond, '_TEST_covariates_02_REGRESSED_', fc, '.png')))
          ### plot PCA of REGRESSED data colored by cur.cond
          png(filename = paste0(test.dir, '/PCA_vst_REGRESSED.', fc, '_col.', cur.cond, '.png'), width = 1100, height = 1000)
          library(ggfortify)
          print(ggplot2::autoplot(prcomp(t(tmp.mat)), data = tmp.annot, colour = cur.cond, size = 3))
          dev.off()
        }
      }
    }
      
      
    ## Remove temporary VST object
    # rm(DE2obj.norm)
    
    ## Performing the DE test ====
    de.try <- try(htg.de.wald <- DESeq2::DESeq(DE2obj))
    if (is(de.try, class2 = 'try-error')) {
      message('DEA failed (see DESeq2 error  above) !')
      file.rename(from = de.dir, to = sub(pattern = '/$', replacement = '_FAILED', x = de.dir))
      next()
    }
    ## Saving the Wald test DESeq object
    # if(save.wald) saveRDS(object = htg.de.wald, file = paste0(factor.dir, '/', cur.cond, '_wald.RDS'), compress = 'bzip2')
    if(save.wald) saveRDS(object = htg.de.wald, file = paste0(de.dir, '/', cur.cond, '_wald.RDS'), compress = 'bzip2')
    
    ## Creating output dir
    mycoef <- paste0(cur.cond, ' : ', paste(cur.condA, collapse='+'), ' -vs- ', paste(cur.condB, collapse='+'))
    message(mycoef)
    message(cur.name)

    # de.dir <- paste(c(factor.dir, cur.name), collapse = '/')
    # dir.create(path = de.dir, recursive = TRUE)

    ## Getting results table for current contrast ====
    mycontrast <- list(paste0(cur.cond, cur.condA), paste0(cur.cond, cur.condB))
    # mycontrast <- c(cur.cond, cur.condA, cur.condB)
    
    # mycontrast <- sapply(all.combz[[mycomb]], function(x) { paste0(cur.cond, x)}, simplify = FALSE)
    # mycontrast <- list(cur.condA, cur.condB)
    # myvalues <- c(rep(1/length(cur.condA), length(cur.condA)), rep(-1/length(cur.condB), length(cur.condB)))
    myvalues <- c(1/length(cur.condA), -1/length(cur.condB))
    DEres <- if (ihw) DESeq2::results(htg.de.wald, contrast = mycontrast, listValues = myvalues, independentFiltering = TRUE, alpha = adjp.max, pAdjustMethod = "BH", parallel = TRUE, BPPARAM = BPPARAM, filterFun = IHW::ihw) else DESeq2::results(htg.de.wald, contrast = mycontrast, listValues = myvalues, independentFiltering = TRUE, alpha = adjp.max, pAdjustMethod = "BH", parallel = TRUE, BPPARAM = BPPARAM)
    ## Saving the test results object
    saveRDS(object = DEres, file = paste0(de.dir, '/', cur.name, '_results.RDS'), compress = 'bzip2')
    
    ## Shrinking l2fc ====
    if (lfcShrink) {
      suppressMessages(DEres <- DESeq2::lfcShrink(htg.de.wald, contrast = mycontrast, type = 'ashr', res = DEres))
      ## Saving the test reults object
      saveRDS(object = DEres, file = paste0(de.dir, '/', cur.name, '_results_lfcShrink.RDS'), compress = 'bzip2')
    }
    # rm(htg.de.wald)
    
    ## Histogram of P-values ====
    png(paste0(de.dir, '/', cur.name, '_phist.png'), width = 2048, height = 768)
    par <- par(mfrow = c(1, 2))
    hist(DEres$pvalue, col = "lightblue", main = paste0("Histogram of raw P-values (DESeq2)\n", mycoef), breaks = 100, xlim = c(0,1), xlab = "P-value")
    hist(DEres$padj, col = "lightblue", main = paste0("Histogram of BH-adjusted P-values (DESeq2)\n", mycoef), breaks = 100, xlim = c(0,1), xlab = "P-value")
    abline(v = adjp.max, col = 2, lty = 2)
    dev.off()
    
    ## MAplot ====
    png(paste0(de.dir, '/', cur.name, '_MA.png'), width = 1024, height = 768)
    DESeq2::plotMA(DEres, alpha = adjp.max, main = paste0("M-A Plot\n", mycoef), cex = 1)
    dev.off()
    
    ## Volcano plot ====
    deg.idx <- DEres$padj <= adjp.max & abs(DEres$log2FoldChange) >= lfc.min
    png(paste0(de.dir, '/', cur.name, '_volcano.png'), width = 1024, height = 768)
    plot(x = DEres$log2FoldChange, y = -log10(DEres$padj), xlab = "log2(Fold-Change)", ylab = "-log10(adjusted P-value)", col = ifelse(deg.idx, "red", "black"), main = paste0("Volcano plot\n", mycoef), pch = 20)
    grid()
    abline(h = -log10(adjp.max), lty = 2, col = 4)
    abline(v = lfc.min * c(-1, 1), lty = 2, col = 4)
    dev.off()
    
    if (length(cur.covars) > 0) norm.mat <- bc.norm.mat
    
    ## Computing per-class metrics ====
    mc.samp.idx <- lapply(levels(SummarizedExperiment::colData(DE2obj)[[cur.cond]]), function(mc.lev) {
      which(SummarizedExperiment::colData(DE2obj)[[cur.cond]] == mc.lev)
    })
    message(table(SummarizedExperiment::colData(DE2obj)[[cur.cond]], useNA = 'always'))
    names(mc.samp.idx) <- levels(SummarizedExperiment::colData(DE2obj)[[cur.cond]])
    mc.metrics <- lapply(names(mc.samp.idx), function(mc.lev) {
      mydf <- data.frame(N = length(mc.samp.idx[[mc.lev]])
                         # , Min = matrixStats::rowMins(norm.mat[,mc.samp.idx[[mc.lev]], drop = FALSE], na.rm = TRUE)
                         , Min = matrixStats::rowMins(x = norm.mat, cols = mc.samp.idx[[mc.lev]], na.rm = TRUE)
                         # , Max = matrixStats::rowMaxs(norm.mat[,mc.samp.idx[[mc.lev]]], drop = FALSE, na.rm = TRUE)
                         , Max = matrixStats::rowMaxs(x = norm.mat, cols = mc.samp.idx[[mc.lev]], drop = FALSE, na.rm = TRUE)
                         , Mean = matrixStats::rowMeans2(x = norm.mat, cols = mc.samp.idx[[mc.lev]], drop = FALSE, na.rm = TRUE)
                         # , Median = matrixStats::rowMedians(x = norm.mat, cols = mc.samp.idx[[mc.lev]], na.rm = TRUE)
                         , Median = matrixStats::rowMedians(x = norm.mat, cols = mc.samp.idx[[mc.lev]], na.rm = TRUE)
                         , matrixStats::rowQuantiles(x = norm.mat, cols = mc.samp.idx[[mc.lev]], probs = c(.25, .75), na.rm = TRUE))
      colnames(mydf) <- paste0(mc.lev, '.', c(colnames(mydf)[1:(ncol(mydf)-2)], 'Q25', 'Q75'))
      return(mydf)
    })
    
    
    ## Output table ====
    ### Bind the feature names as 1st column, then ID conversions, then DEresults, then metrics
    DEres.df <- cbind(Feature = rownames(DEres), sconv[rownames(DEres),], as.data.frame(DEres), Reduce(f = cbind, x = mc.metrics))
    
    sig.word <- paste0('Sig_@adjp', adjp.max, '_lfc', lfc.min)
    DEres.df[[sig.word]] <- 0
    DEres.df[[sig.word]][deg.idx] <- 1
    
    ## Sorting a temp df for export
    DEres.df.out <- DEres.df[order(DEres.df$padj, abs(DEres.df$log2FoldChange), decreasing = c(FALSE, TRUE)),]
    write.table(DEres.df.out, file = paste0(de.dir, '/', cur.name, '_results.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
    
    # sig.genes <- as.character(DEres.df$Feature[DEres.df[sig.word] == 1])
    # sig.genes <- as.character(DEres.df$SYMBOL[DEres.df[sig.word] == 1])
    sig.feats <- which(DEres.df[sig.word] == 1) ## Line positions
    
    # ## Setting SYMBOLS as norm.mat ids (to fit with sig.genes)
    # rownames(norm.mat) <- sconv$SYMBOL[rownames(norm.mat)]
    
    ## Draw feature boxplots ====
    # if(boxplots & length(sig.genes) > 0) {
    if(boxplots & length(sig.feats) > 0) {
      boxdir <- paste0(de.dir, '/boxplots')
      dir.create(path = boxdir, recursive = TRUE)
      # for (g in sig.genes[1:(min(length(sig.genes), or.top.max))]) {
      for (g in sig.feats[1:(min(length(sig.feats), or.top.max))]) {
        # message('Feature : ', DEres.df$Feature[g])
        # message('\t', DEres.df$Feature[g] %in% rownames(norm.mat))
        png(filename = paste0(boxdir, '/', DEres.df$SYMBOL[g], '_norm.exp_boxplot.png'), width = 800, height = 600)
        sig.split <- split(norm.mat[g, ], droplevels(SummarizedExperiment::colData(DE2obj)[[cur.cond]]))
        boxplot2l(x = sig.split, col = seq_along(sig.split), main = paste0(DEres.df$SYMBOL[g], ' normalized expression VS ', cur.cond, '\nDESeq2 : l2FC = ', round(DEres.df[g, 'log2FoldChange'], digits = 3), ' ; adjP = ', format(DEres.df[g, 'padj'], scientific = TRUE, digits = 3)), xlab = '', ylab = 'Normalized expression', vertical = TRUE, cex = 2)
        
        dev.off()
      }
    }
    
    ## Setting a color palette for the heatmaps
    myRamp <- circlize::colorRamp2(c(-2, 0, 2), heatmap.palette)
      
    # if (length(sig.genes) > enr.min.genes) {
    if (length(sig.feats) > enr.min.genes) {
        
      cur.annot <- as.data.frame(SummarizedExperiment::colData(DE2obj)[,c(cur.cond, cur.covars), drop = FALSE])
      
      ## Heatmap
      ## data to plot
      # plotDat <- norm.mat[rownames(norm.mat) %in% sig.genes,]
      plotDat <- norm.mat[sig.feats,]
      z.mat <- (plotDat - rowMeans(plotDat)) / matrixStats::rowSds(plotDat)
      # Creating sample annotation
      set.seed(1)
      ha1 = ComplexHeatmap::HeatmapAnnotation(df = cur.annot[,c(cur.cond, cur.covars), drop = FALSE])
      
      ## Looping through requested clustering methods
      for (sdm in samples.dist.method) {
        for (shm in samples.hclust.method) {
          for (gdm in genes.dist.method) {
            for (ghm in genes.hclust.method) {
              ## Clustering samples
              hc.s <- hclust(amap::Dist(x = t(plotDat), method = sdm), method = shm)
              ## Clustering genes
              hc.g <- hclust(amap::Dist(x = plotDat, method = gdm), method = ghm)
              ## Compute heatmap ====
              set.seed(my.seed)
              myHM <- suppressMessages(ComplexHeatmap::Heatmap(z.mat, name = "Normalized counts"
                                                               # use my custom color palette
                                                               , col = myRamp
                                                               # do not show gene names
                                                               , show_row_name = TRUE
                                                               # do not clusterize samples
                                                               , cluster_columns = hc.s
                                                               , cluster_rows = hc.g
                                                               # add a nice grey border to cells
                                                               , rect_gp = grid::gpar(col = "darkgrey", lwd=0.5)
                                                               # add sample annotation
                                                               , top_annotation = ha1
                                                               , use_raster = TRUE
                                                               , raster_device = 'png'
                                                               ))
              ## Draw heatmap ====
              png(paste0(de.dir, '/', cur.name, '_sig.', nrow(z.mat), 'x', ncol(z.mat), '_', paste(c(gdm, ghm, sdm, shm), collapse = "_"), '.heatmap.png'), width = min(ncol(z.mat) * 15, 2000) + 600, height = min(length(sig.feats) * 10, 5000) + 300)
              ComplexHeatmap::draw(myHM)
              dev.off()
            }
          }
        }
      }
    }
    
  
    ### Shorter heatmap if more sig genes than requested "topN" ====
    # if(length(sig.genes) > or.top.max) {
    if(length(sig.feats) > or.top.max) {
        
      ## Heatmap
      ## data to plot
      # sig.genes <- as.character(DEres.df$Feature[DEres.df[[sig.word]] == 1][1:or.top.max])
      sig.feats <- sig.feats[1:or.top.max]
      # plotDat <- norm.mat[rownames(norm.mat) %in% sig.genes,]
      plotDat <- norm.mat[sig.feats,]
      z.mat <- (plotDat - rowMeans(plotDat)) / matrixStats::rowSds(plotDat)
      # Creating sample annotation
      set.seed(1)
      ha1 = ComplexHeatmap::HeatmapAnnotation(df = cur.annot[,c(cur.cond, cur.covars), drop = FALSE])
      ## Clustering samples
      for (sdm in samples.dist.method) {
        for (shm in samples.hclust.method) {
          for (gdm in genes.dist.method) {
            for (ghm in genes.hclust.method) {
              hc.s <- hclust(amap::Dist(x = t(plotDat), method = sdm), method = shm)
              ## Clustering genes
              hc.g <- hclust(amap::Dist(x = plotDat, method = gdm), method = ghm)
              ## Computing heatmap
              myHM <- suppressMessages(ComplexHeatmap::Heatmap(z.mat, name = "Normalized counts",
                                                             col = myRamp,
                                                             show_row_name = TRUE,
                                                             cluster_columns = hc.s,
                                                             cluster_rows = hc.g,
                                                             rect_gp = grid::gpar(col = "darkgrey", lwd=0.5),
                                                             top_annotation = ha1))
              ## Draw
              png(paste0(de.dir, '/', cur.name, '_sig.TOP', nrow(z.mat), 'x', ncol(z.mat), '_', paste(c(gdm, ghm, sdm, shm), collapse = "_"), '.heatmap.png'), width = min(ncol(z.mat) * 15, 2000) + 600, height = min(length(sig.feats) * 10, 5000) + 300)
              ComplexHeatmap::draw(myHM)
              dev.off()
            }
          }
        }
      }
    }
  
    ### Functional enrichment ====
    if (any(msigdb.do, kegg.do, do.do, go.do, wp.do, reactome.do, mesh.do) & !is.null(species)) {
      
      if(length(sig.feats) < enr.min.genes) {
        message(paste0('Less significant genes (', length(sig.feats), ') at the defined thresholds than the minimum expected (', enr.min.genes, ') : ORA analyses are deactivated.'))
        msigdb.do[2] <- kegg.do[2] <- do.do[2] <- go.do[2] <- wp.do[2] <- reactome.do[2] <- mesh.do[2] <- FALSE
      }
      
      enr.inputs <- table2enr(deseq2.res.data = DEres.df, species = species, geneid.colname = 'Feature', geneid.type = feature.type, value.colname = 'log2FoldChange', topN.max = or.top.max, topN.order.colname = 'padj', topN.order.decreasing = FALSE, topN.cutoff = enrp.max, topN.keep.operator = '<')
        
      ## MSIGDB
      if (any(msigdb.do)) {
        message('\t[MSigDb]')
        msigdb.collec <- as.data.frame(msigdbr::msigdbr_collections())
        ## Removing redundancy when requested
        if(non.redundant) {
          if (any(go.do)) msigdb.collec <- msigdb.collec[-c(grep(pattern = 'GO\\:', x = msigdb.collec$gs_subcat)),]
          if (any(wp.do)) msigdb.collec <- msigdb.collec[-c(grep(pattern = 'WIKIPATHWAYS', x = msigdb.collec$gs_subcat)),]
          if (any(reactome.do)) msigdb.collec <- msigdb.collec[-c(grep(pattern = 'REACTOME', x = msigdb.collec$gs_subcat)),]
          if (any(kegg.do)) msigdb.collec <- msigdb.collec[-c(grep(pattern = 'KEGG', x = msigdb.collec$gs_subcat)),]
        }
        for (mx in seq_len(nrow(msigdb.collec))) {
          msc <- if(msigdb.collec[mx,2] == '') msigdb.collec[mx,1] else paste(c(msigdb.collec[mx,1], msigdb.collec[mx,2]), collapse = '_')
          ## Import the TERM2GENE object corresponding to the desired category/subcategory combo
          my.t2g <- msigdb_to_t2g(species = species, category = msigdb.collec[mx,1], subcategory = msigdb.collec[mx,2])
          ### GSEA
          if (msigdb.do[1]) {
            my.gsea.res <- try(gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = 'clusterProfiler::GSEA', t2g = my.t2g, t2g.name = msc, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
            ## Generate plots / outputs
            if (!is(my.gsea.res, class2 = 'try-error')) gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = cur.name)
          }
          ### ORA
          if(msigdb.do[2]) {
            my.ora.res <- try(ora.run(gene = enr.inputs$ora.genevec, universe = unname(enr.inputs$gene2Symbol), species = species, func.name = 'clusterProfiler::enricher', t2g = my.t2g, t2g.name = msc, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
            ## Generate plots / outputs
            if (!is(my.ora.res, class2 = 'try-error')) ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = cur.name, geneList = enr.inputs$gsea.genevec)
          }
        }
      }
      
      ## GO (gene ontology)
      if (any(go.do)) {
        message('\t[GO]')
        ### GSEA
        if(go.do[1]) {
          func.name <- 'clusterProfiler::gseGO'
          for (x in c('BP', 'CC', 'MF')) {
            my.org <- paste0(msigdbr2org(species), '.db')
            library(my.org, character.only = TRUE)
            my.gsea.res <- try(gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enrp.max, minGSSize = enr.min.genes, OrgDb = get(my.org), ont = x))
            if (!is(my.gsea.res, class2 = 'try-error')) {
              my.gsea.res@setType <- paste(c(my.gsea.res@setType, x), collapse = '_')
              gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = cur.name)
              ## Simplify
              if(nrow(my.gsea.res) > 1) {
                my.gsea.res@setType <- x
                my.gsea.res <- enrichplot::pairwise_termsim(my.gsea.res)
                my.gsea.res.simp <- clusterProfiler::simplify(my.gsea.res, cutoff = 0.7, by = "p.adjust", select_fun = min)
                if(nrow(my.gsea.res.simp) < nrow(my.gsea.res)) {
                  my.gsea.res.simp@setType <- paste(c(func.name, x, 'simplified'), collapse = '_')
                  gsea.output(gseaResult = my.gsea.res.simp, out.dir = de.dir, comp.name = cur.name)
                }
              }
            }
          }
        }
        ### ORA
        if(go.do[2]) {
          func.name <- 'clusterProfiler::enrichGO'
          for (x in c('BP', 'CC', 'MF')) {
            my.org <- paste0(msigdbr2org(species), '.db')
            library(my.org, character.only = TRUE)
            my.ora.res <- ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes, OrgDb = get(my.org), ont = x)
            if (!is(my.ora.res, class2 = 'try-error')) {
              my.ora.res@ontology <- paste(c(my.ora.res@ontology, x), collapse = '_')
              ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = cur.name, geneList = enr.inputs$gsea.genevec)
              ## Simplify
              if(nrow(my.ora.res) > 1) {
                my.ora.res@ontology <- x
                my.ora.res <- enrichplot::pairwise_termsim(my.ora.res)
                my.ora.res.simp <- clusterProfiler::simplify(my.ora.res, cutoff = 0.7, by = "p.adjust", select_fun = min)
                if(nrow(my.ora.res.simp) < nrow(my.ora.res)) {
                  my.ora.res.simp@ontology <- paste(c(func.name, x, 'simplified'), collapse = '_')
                  ora.output(enrichResult = my.ora.res.simp, out.dir = de.dir, comp.name = cur.name, geneList = enr.inputs$gsea.genevec)
                }
              }
            }
          }
        }
      }
      
      ## DO (disease ontology)
      if (any(do.do)) {
        message('\t[DO]')
        ### GSEA
        if (do.do[1]) {
          for (x in c('DOSE::gseDO', 'DOSE::gseNCG', 'DOSE::gseDGN')) {
            my.gsea.res <- try(gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
            if (!is(my.gsea.res, class2 = 'try-error')) gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = cur.name)
          }
          ### ORA
          if(do.do[2]) {
            for (x in c('DOSE::enrichDO', 'DOSE::enrichNCG', 'DOSE::enrichDGN')) {
              my.ora.res <- try(ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
              if (!is(my.ora.res, class2 = 'try-error')) ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = cur.name, geneList = enr.inputs$gsea.genevec)
            }
          }
        }
      }
      
      ## KEGG/MKEGG
      ### NOTE1 : It's the same way to call the 'gsea.run' / 'ora.run' as it is for 'DO', 'NCG' or 'DGN', but here it's compatible with many more species than homo sapiens.
      ### NOTE2 : for this case, additional KEGG pathway plots will be generated.
      ### NOTE3 : for this case, an internet connexion is required to query the KEGG website.
      if (any(kegg.do)) {
        message('\t[KEGG]')
        ### GSEA
        if (kegg.do[1]) {
          for (x in c('clusterProfiler::gseKEGG', 'clusterProfiler::gseMKEGG')) {
            my.gsea.res <- try(gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
            if (!is(my.gsea.res, class2 = 'try-error')) gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = cur.name)
          }
          ### ORA
          if(kegg.do[2]) {
            for (x in c('clusterProfiler::enrichKEGG', 'clusterProfiler::enrichMKEGG')) {
              my.ora.res <- try(ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
              if (!is(my.ora.res, class2 = 'try-error')) ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = cur.name, geneList = enr.inputs$gsea.genevec)
            }
          }
        }
      }
      
      ## WP (wikipathways)
      if(any(wp.do)) { 
        message('\t[WP]')
        if(wp.do[1]) {
          ### GSEA
          func.name <- 'clusterProfiler::gseWP'
          my.gsea.res <- try(gsea.run(geneList = enr.inputs$gsea.genevec, organism = species, func.name = func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
          if (!is(my.gsea.res, class2 = 'try-error')) gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = cur.name)
        }
        ### ORA
        if(wp.do[2]) {
          func.name <- 'clusterProfiler::enrichWP'
          my.ora.res <- ora.run(gene = enr.inputs$ora.genevec, organism = species, func.name = func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes)
          if (!is(my.ora.res, class2 = 'try-error')) ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = cur.name, geneList = enr.inputs$gsea.genevec)
        }
      }
      
      ## REACTOME
      if (any(reactome.do)) {
        message('\t[Reactome]')
        org.name <- paste0(msigdbr2org(species = species), '.db')
        library(org.name, character.only = TRUE)
        reactome.org <- tolower(convert_species_name(OrgDb = get(org.name)))
        if(reactome.do[1]) {
          ### GSEA
          func.name <- 'ReactomePA::gsePathway'
          my.gsea.res <- gsea.run(geneList = enr.inputs$gsea.genevec, organism = reactome.org, func.name = func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enrp.max, minGSSize = enr.min.genes)
          my.gsea.res@setType <- paste0(func.name, '_Reactome')
          gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = cur.name)
        }
        if(reactome.do[2]) {
          ### ORA
          func.name <- 'ReactomePA::enrichPathway'
          my.ora.res <- ora.run(gene = enr.inputs$ora.genevec, organism = reactome.org, func.name = func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes)
          my.ora.res@ontology <- paste0(func.name, '_Reactome')
          ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = cur.name, geneList = enr.inputs$gsea.genevec)
        }
      }
      
      ## MESH (WARNING : MEMORY OGRE AND SLOW !! Big DBs, 3 sources, 16 categories ! 64 GB of RAM required for most bases !
      ### Requires additional parameters :
      ### . 'MeSHDb' : character ; name of a MeSH [NO : AUTO FROM SPECIES NAME]
      ### . 'database' : character ; MeSH source type (can be 'gendoo' = text-mining, 'gene2pubmed' = manual curation by NCBI team, 'RBBH' = sequence homology with BLASTP search @ E-value < 1E-50)
      ### . 'category' : character ; name of a MeSH category sub-db (namely 'A', 'B', 'C', 'D', 'G').
      ### NOTE : see https://yulab-smu.top/biomedical-knowledge-mining-book/meshes-semantic-similarity.html
      if (any(mesh.do)) {
        message('\t[MeSHDb]')
        ### List of requested MeSH DBs
        mesh.dbs <- c('gendoo', 'gene2pubmed', 'RBBH') ## 'RBBH' is not available for Homo sapiens.
        ### List of requested MeSH categories
        mesh.categories <- toupper(letters[-c(15:21,23:25)]) ## More categories are available, but some do not seem to work with Homo sapiens for some of the DBs.
        ### Building the MeSH package name corresponding to the current species
        mesh.sp <- paste0(c('MeSH.', substr(unlist(strsplit(species, ' ')), c(1, 1), c(1,2)), '.eg.db'), collapse = '')
        ### Checking which MeSH DBs are available for the current species.
        mesh.dbs <- MeSHDbi::listDatabases(eval(parse(text = paste0(mesh.sp, '::', mesh.sp))))[,1]
        
        ### ORA
        #### WARNING !! the 'gene2pubmed' requires a lot of RAM (~12 GB) !!
        if (mesh.do[2]) {
          mesh.func.name <- 'meshes::enrichMeSH'
          for (y in mesh.dbs) {
            for (x in mesh.categories) {
              message(paste0(y, ' ', x))
              if (y %in% mesh.dbs) {
                my.ora.res <- try(ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = mesh.func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes, database = y, category = x))
                if(!is(my.ora.res, class2 = 'try-error')) {
                  ## Little hack specific to MeSH results (as I was not able to get the value of extra parameters 'database' and 'category' from within the 'gsea.run()' function)
                  my.ora.res@ontology <- paste(c(mesh.func.name, y, x), collapse = '_')
                  ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = cur.name, geneList = enr.inputs$gsea.genevec)
                }
              } else message(paste0("Unsupported MeSH database '", y, "'. Expecting one of : '", paste(mesh.dbs, collapse = "', '"), "'."))
            }
          }
        }
        
        ### GSEA
        #### WARNING !! Needs too much memory for a laptop (probably over 64 GB of RAM, easily...). SO, not recommended out of flamingo.
        if(mesh.do[1]) {
          mesh.func.name <- 'meshes::gseMeSH'
          for (y in mesh.dbs) {
            for (x in mesh.categories) {
              message(paste0(y, ' ', x))
              if (y %in% mesh.dbs) {
                my.gsea.res <- try(gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = mesh.func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enrp.max, minGSSize = enr.min.genes, database = y, category = x), silent = TRUE)
                if (!is(my.gsea.res, class2 = 'try-error')) {
                  ## Little hack specific to MeSH results (as I was not able to get the value of extra parameters 'database' and 'category' from within the 'gsea.run()' function)
                  my.gsea.res@setType <- paste(c(mesh.func.name, y, x), collapse = '_')
                  gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = cur.name)
                }
              } else message(paste0("Unsupported MeSH database '", y, "'. Expecting one of : '", paste(mesh.dbs, collapse = "', '"), "'."))
            }
          }
        }
      }
    }
  }
}
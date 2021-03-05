################
#### HEADER ####
################

## A Wrapper to automate GSEA/ORA functional annotation and analysis, based on clusterProfiler
## R [4.0.4]
## PACKAGES :
### REQUIRED :
#### CRAN : 'tidyverse' [1.3.0], 'ggnewscale' [0.4.5], 'vroom' [1.4.0]
#### Bioconductor : 'clusterProfiler' [3.18.1] (will import other required packages as dependencies, like 'DOSE' and 'enrichplot'), 'org.Xx.eg.db' (species-specific, like 'org.Hs.eg.db' [3.12.0])
### SUGGESTED :
#### Bioconductor : 'msigdbr' [7.2.1] (to assess MSigDB databases), 'meshes' [1.16.0] (to assess MeSH databases), 'MeSH.Xxx.eg.db' (species-specific, like 'MeSH.Hsa.eg.db' [1.15.0], required to query MeSH databases), 'pathview' [1.30.1] (to plot KEGG pathways if 'clusterProfiler::enrichKEGG' or 'clusterProfiler::gseKEGG' functions are used)



###################
#### FUNCTIONS ####
###################

### FUNCTION TO CONVERT AN msigdbr SPECIES NAME TO AN AnnotationForge ORGANISM PACKAGE ROOTNAME (ie, the package name without the '.db' suffix)
### . 'species' : character ; a species name, as in the 'species_name' column of msigdbr::msigdbr_species()
### NOTE1 : This converter should work with all available species in msigdbr, with the exception of saccharomyces cerevisiae (which does not have 'eg' in its organism name)
### NOTE 2 : See MSigDB collections description here : http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
msigdbr2org <- function(species = NULL) {
  return(paste0('org.', paste0(stringr::str_sub(unlist(base::strsplit(species, ' ')), 1, 1), collapse = ''), '.eg'))
}

### FUNCTION TO LOAD AN MSigDB BANK AND CREATE THE term2gene TABLE
### . 'species' : character ; a species name, as in the 'species_name' column of msigdbr::msigdbr_species()
### . 'category' : character ; an MSigDB category, as in the 'gs_cat' column of msigdbr::msigdbr_collections()
### . 'subcategory' : character ; an MSigDB category, as in the 'gs_subcat' column of msigdbr::msigdbr_collections()
msigdb_to_t2g <- function(species = 'Homo sapiens', category = NULL, subcategory = NULL) {
  ## Checking species
  valid.species <- msigdbr::msigdbr_species()
  if (!species %in% valid.species$species_name) stop(paste0("Unsupported species. Expecting any of '", paste(valid.species$species_name, collapse = "', '"), "'. See ?msigdbr::msigdbr"))
  ## Getting the bank data
  message(paste0('Loading MSigDB data for ', category, ' ', subcategory))
  msigdb <- msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)
  ## Generating the output table
  t2g <- data.frame(term = as.factor(msigdb$gs_name), gene = as.character(msigdb$entrez_gene), stringsAsFactors = FALSE)
  return(t2g)
}

## Function to create input objects for gsea.run() and ora.run() from an output of our rna-salmon-deseq2 pipeline
### . 'deseq2.res.file' : character ; output table from the rna-salmon-deseq2 pipeline (the *_complete.tsv table)
### . 'species' : character ; species name (Homo sapiens, Mus musculus, etc ...)
### . 'geneid.colname' : character ; column name in 'deseq2.res.file' input for the gene identifier.
### . 'geneid.type' : character ; type of gene identifier. Should be one value of the output of clusterProfiler::idType("org.Xx.eg.db") (with 'Xx' corresponding as the species of interest). Usually 'SYMBOL' or 'ENSEMBL'.
### . 'value.colname' : character ; column name in 'deseq2.res.file' input for the (numeric/integer) values to use in GSEA (ie, logFoldChange).
### . 'topN.max' : integer ; maximum number of significant genes to keep for ORA (may be inferior if fewer genes were found as differentially expressed)
### . 'topN.order.colname' : character ; column name in 'deseq2.res.file' input to use to order the table before selecting the 'topN' genes for ORA.
### . 'topN.order.decreasing' : logic ; use decreasing order for 'topN.order.colname'.
### . 'topN.cutoff' : numeric ; cutoff to use to select 'topN' genes on 'topN.order.colname' when there are less significant genes than 'topN'.
### . 'topN.keep.operator' : character ; operator (given as a character) to keep genes when evaluating their 'topN.order.colname' value to 'topN.cutoff' (usually one of '<', '<=', '>', '>=').
### . '...' : any parameter to read.table() to handle the input file (usually 'sep = "\t", header = TRUE, as.is = TRUE')
pipe2enr <- function(deseq2.res.file = NULL, species = "Homo sapiens", geneid.colname = 'Gene_Name', geneid.type = 'SYMBOL', value.colname = 'stat_change', topN.max = 100, topN.order.colname = 'Adjusted_PValue', topN.order.decreasing = FALSE, topN.cutoff = 5E-02, topN.keep.operator = '<', ...) {
  ## Checks
  ## Loading the DE table
  deres <- read.table(file = deseq2.res.file, ...)
  ## Forcing gene ID to be a character
  deres[[geneid.colname]] <- as.character(deres[[geneid.colname]])
  ## Converting non-Entrez IDs to Entrez IDs
  if (geneid.type != 'ENTREZID') {
    ## Building the gene ID conversion vector
    gconv <- as.list(clusterProfiler::bitr(deres[[geneid.colname]], fromType=geneid.type, toType='ENTREZID', OrgDb=paste0(msigdbr2org(species), '.db')))
    names(gconv$ENTREZID) <- gconv[[geneid.type]]
    names(gconv[[geneid.type]]) <- gconv$ENTREZID
    deres$ENTREZID <- gconv$ENTREZID[as.character(deres[[geneid.colname]])]
  } else {
    gconv <- list(ENTREZID = setNames(deres$ENTREZID, deres$ENTREZID))
    deres$ENTREZID <- deres[[geneid.colname]]
  }
  ## Removing entries without EntrezID
  na.entrez <- is.na(deres$ENTREZID)
  deres <- deres[!na.entrez,]
  message(paste0('Removed ', length(which(na.entrez)), ' lines without EntrezID value.'))
  ## Removing entries with duplicated EntrezID
  dup.entrez.names <- unique(deres$ENTREZID[duplicated(deres$ENTREZID)])
  dup.entrez <- deres$ENTREZID %in% dup.entrez.names
  deres <- deres[!dup.entrez,]
  message(paste0('Removed ', length(which(dup.entrez)), ' lines (', length(dup.entrez.names), ' genes) with replicated EntrezID value.'))
  ## Formatted input
  ## GSEA
  gsea.genevec <- sort(setNames(deres[[value.colname]], deres$ENTREZID), decreasing = TRUE)
  ## ORA
  deres2 <- deres[!is.na(deres[[topN.order.colname]]),]
  deres2 <- deres2[order(deres2[[topN.order.colname]], decreasing = topN.order.decreasing),]
  top.keep <- min(length(which(eval(parse(text = paste0('deres2[["', topN.order.colname, '"]] ', topN.keep.operator, ' topN.cutoff'))))), topN.max)
  ora.genevec <- if(length(top.keep) > 0) setNames(deres2$ENTREZID[1:top.keep], deres2[[geneid.colname]][1:top.keep]) else c()
  
  return(list(gsea.genevec = gsea.genevec,
              ora.genevec = ora.genevec,
              gene2Symbol = gconv$ENTREZID))
}

### FUNCTION TO PERFORM GSEA
### . 'geneList' : numeric vector : a named vector of decreasing values, with EntrezIDs as names
### . 'species' : character ; a species name, as in the 'species_name' column of msigdbr::msigdbr_species()
### . 'category' : character ; an MSigDB category, as in the 'gs_cat' column of msigdbr::msigdbr_collections()
### . 'subcategory' : character ; an MSigDB category, as in the 'gs_subcat' column of msigdbr::msigdbr_collections()
### . 'gene2Symbol' : character vector ; a named vector of EntrezIDs, with corresponding Symbol or EnsemblID as names
### . 'seed' : integer ; the random number generation seed
### . 'pvalueCutoff' : numeric ; minimum FDR-adjusted p-value for signficantly enriched terms (see clusterProfiler::GSEA)
### . 'minGSSize' : integer ; minimal size of each geneSet for analyzing (see clusterProfiler::GSEA)
### . '...' : any other parameter to pass to 'func.name'
### NOTE : For input ('geneList'), ALWAYS use all available genes (ie, not limited to significant differentialy expressed genes), as intended for GSEA analysis. NEVER use a selection of genes (results will be corrupt). Please also be careful that you may obtain significant GSEA results from a list of values corresponding to NO significant genes (as GSEA only needs the order of these genes). In this case you may talk about "tendencies" of enrichment in your differential expression results.
gsea.run <- function(geneList = NULL, func.name = 'clusterProfiler::GSEA', species = 'Homo sapiens', t2g = NULL, t2g.name = NULL, gene2Symbol = NULL, seed = 1337, pvalueCutoff = 5E-02, minGSSize = 10, ...) {

  ## Checks
  if (any(!is.numeric(geneList))) stop("'geneList' should be a decreasingly sorted numeric vector, named with EntrezIDs.")
  if (length(names(geneList)) == 0) stop("'geneList' should be a decreasingly sorted numeric vector, named with EntrezIDs.")
  if (!all(geneList == sort(geneList, decreasing = TRUE))) stop("'geneList' should be a decreasingly sorted numeric vector, named with EntrezIDs.")
  if (!is.character(func.name)) stop("'func.name' should be a character.")
  if (!is.character(species)) stop("'species' should be a character.")
  valid.species <- msigdbr::msigdbr_species()$species_name
  if (!(species %in% valid.species)) stop(paste0("Unsupported 'species'. Expecting one of : '", paste(valid.species, collapse = "', '"), '.'))
  if (!is.null(t2g) & is.null(t2g.name)) stop("'t2g' provided but no 't2g.name'.")
  if (!is.null(t2g.name) & is.null(t2g)) stop("'t2g.name' provided but no 't2g'.")
  if (!is.null(gene2Symbol)) {
    if (any(!is.character(gene2Symbol))) stop("'gene2Symbol' shoud be a character vector, named with EntrezIDs.")
    if (length(names(gene2Symbol)) == 0) stop("'gene2Symbol' shoud be a character vector, named with EntrezIDs.")
  }
  if (!is.numeric(seed)) stop("'seed' should be an integer.")
  if (!is.numeric(pvalueCutoff)) stop("'pvalueCutoff' should be a positive, < 1, numeric")
  if (!(pvalueCutoff >= 0)) stop("'pvalueCutoff' should be a positive, < 1, numeric")
  if (!(pvalueCutoff <= 1)) stop("'pvalueCutoff' should be a positive, < 1, numeric")
  if (!is.numeric(minGSSize)) stop("'minGSSize' should be an integer.")
  if (!(minGSSize > 0)) stop("'minGSSize' should a non-null positive integer.")
  if (minGSSize < 10) warning("'minGSSize' < 10 : expect no result !")
  
  ## Loading the requested function
  func.split <- unlist(strsplit(func.name, '::'))
  gse.function <- base::get(func.split[2], envir = loadNamespace(func.split[1]))

  if (func.name == 'clusterProfiler::GSEA') {
    if (is.null(t2g)) stop("When calling this function with func.name = 'clusterProfiler::GSEA', a value is required for 't2g'")
    if (is.null(t2g.name)) stop("When calling this function with func.name = 'clusterProfiler::enricher', a value is required for 't2g.name'")
    
    ## Functions requiring a TERM2GENE (msigdbr, CellMarkers, ...)
    gsea.res <- gse.function(geneList = geneList, TERM2GENE = t2g, pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, seed = seed, ...)
    gsea.res@setType <- paste(c(func.name, t2g.name), collapse = '_')
  } else if (func.name == 'meshes::gseMeSH') {
    ## MeSH (requires additional 'MeSHDb', 'database' and 'category' parameters)
    mesh.sp <- paste0(c('MeSH.', substr(unlist(strsplit(species, ' ')), c(1, 1), c(1,2)), '.eg.db'), collapse = '')
    gsea.res <- gse.function(geneList = geneList, MeSHDb = mesh.sp, pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, seed = seed, ...)
    # args.all <- as.list(match.call(expand.dots = FALSE))
    gsea.res@setType <- func.name
  } else {
    if ('kegg' %in% tolower(func.name)) {
      ## KEGG / KEGGM (requires a custom species name in 'organism' parameter)
      kegg.sp <- tolower(paste0(substr(unlist(strsplit(species, ' ')), c(1, 1), c(1,2)), collapse = ''))
      gsea.res <- gse.function(geneList = geneList, organism = kegg.sp, pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, seed = seed, ...)
    } else if (species != 'Homo sapiens') stop(paste0("Function '", func.name, "' can only be used with the 'Homo sapiens' species.")) else {
      ## OTHER (generic functions without TERM2GENE or custom parameters, like those in DOSE package)
      gsea.res <- gse.function(geneList = geneList, pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, seed = seed, ...)
    }
    gsea.res@setType <- func.name
  }

  ## Adding some useful metadata
  gsea.res@organism <- species
  if(!is.null(gene2Symbol)) {
    gsea.res@gene2Symbol <- gene2Symbol
    gsea.res@keytype <- 'ENTREZID'
  }
  return(gsea.res)
}

### FUNCTION TO OUTPUT TABLE AND PLOT RESULTS FROM gsea.run()
### . 'gseaResult' : gseaResult object ; an output from the gsea.run() function
### . 'comp.name' : character ; The name of the differential expression comparison which resulted in the 'geneList' input to gsea.run(). Actually, only used in plots' title.
### . 'out.dir' : character ; directory to output plots and tables
### . 'heatplot' : integer of NULL ; if integer perform an enrichplot::dotplot() limited to a number of term corresponding to the given integer value. If null, no plot is performed.
### . 'dotplot' : integer of NULL ; if integer perform an enrichplot::dotplot() limited to a number of term corresponding to the given integer value. If null, no plot is performed.
### . 'barplot' : integer of NULL ; if integer perform an enrichplot::barplot() limited to a number of term corresponding to the given integer value. If null, no plot is performed.
### gsea.plot : logic ; perform per-term enrichment plots like generated by the BROAD GSEA app.
### . 'ridgeplot' : integer of NULL ; if integer perform an enrichplot::ridgeplot() limited to a number of term corresponding to the given integer value. If null, no plot is performed.
### . 'cnetplot' : integer of NULL ; if integer perform an enrichplot::cnetplot() limited to a number of term corresponding to the given integer value. If null, no plot is performed.
### . 'emapplot' : integer of NULL ; if integer perform an enrichplot::emapplot() limited to a number of term corresponding to the given integer value. If null, no plot is performed.
gsea.output <- function(gseaResult = NULL, comp.name = 'TEST', out.dir = getwd(), heatplot = 100, dotplot = 100, barplot = 100, gsea.plot = TRUE, ridgeplot = 100, cnetplot = 10, emapplot = 10) {

  ## Significant terms ?
  gsea.sig.tf <- gseaResult@result$p.adjust < gseaResult@params$pvalueCutoff

  if (any(gsea.sig.tf)) {

    ## Creating GSEA output dir
    gsea.dir <- paste(c(out.dir, paste0('GSEA_adjp.', gseaResult@params$pvalueCutoff), gseaResult@setType), collapse = '/')
    dir.create(path = gsea.dir, recursive = TRUE)

    ## Dumping GSEA results
    saveRDS(gseaResult, file = paste0(gsea.dir, '/GSEA.results.RDS'), compress = 'bzip2')
    
    ## Converting to readable symbols
    gseaResult_readable <- DOSE::setReadable(gseaResult, OrgDb = paste0(msigdbr2org(gseaResult@organism), '.db'))
    ## Writing readable table
    write.table(gseaResult_readable@result, file = paste0(gsea.dir, '/GSEA.results_readable.txt'), quote = FALSE, sep = "\t", row.names = FALSE)

    ## Heatplot (needs readable)
    if (!is.null(heatplot)) {
      hp <- enrichplot::heatplot(gseaResult_readable, showCategory = heatplot, foldChange = gseaResult@geneList)
      png(paste0(gsea.dir, '/GSEA.heatplot.png'), width = (5 * nrow(hp$data)) + 500, height = (min(length(which(gsea.sig.tf)), heatplot) * 12) + 150)
      print(hp)
      dev.off()
    }
    
    ## Dotplot // Barplot
    complot <- c()
    if (!is.null(dotplot)) complot <- c(complot, setNames(dotplot, 'clusterProfiler::dotplot'))
    if (!is.null(barplot)) complot <- c(complot, setNames(barplot, 'enrichplot:::barplot.enrichResult'))
    for (cx in seq_along(complot)) {
      func.split <- unlist(strsplit(names(complot)[cx], ':+'))
      plot.function <- base::get(func.split[2], envir = loadNamespace(func.split[1]))
      gseaBak <- gseaResult
      gseaResult@result$Description <- gsub(pattern = '_', ' ', gseaResult@result$Description, fixed = TRUE)
      png(paste0(gsea.dir, '/GSEA.', func.split[2], '.png'), width = 1000, height = (min(length(which(gsea.sig.tf)), unname(complot[cx])) * 20) + 300)
      print(plot.function(gseaResult, showCategory = unname(complot[cx]), title = paste(c(comp.name, gseaResult@setType, 'GSEA'), collapse = "\n"), split = '.sign') + ggplot2::facet_grid(~ .sign) + ggplot2::scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=100)))
      dev.off()
      gseaResult <- gseaBak
      rm(gseaBak)
    }
    
    if (!is.null(ridgeplot)) {
      ## Ridgeplot
      gseaBak <- gseaResult
      gseaResult@result$Description <- gsub(pattern = '_', ' ', gseaResult@result$Description, fixed = TRUE)
      rip <- enrichplot::ridgeplot(gseaResult, showCategory = ridgeplot)
      png(paste0(gsea.dir, '/GSEA.ridgeplot.png'), width = 1000, height = (min(length(which(gsea.sig.tf)), ridgeplot) * 20) + 300)
      print(rip + ggplot2::ggtitle(paste(c(comp.name, gseaResult@setType, 'GSEA'), collapse = "\n")) + ggplot2::scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=100)))
      dev.off()
      gseaResult <- gseaBak
      rm(gseaBak)
    }
    
    ## GSEA plots
    if (gsea.plot) {
      gseaplot.dir <- paste0(gsea.dir, '/gseaplot_clusterProfiler/')
      dir.create(path = gseaplot.dir, recursive = TRUE)
      png(paste0(gseaplot.dir, '/GSEA.gseaplot_cP_%04d.png'), width = 1024, height = 768)
      for (x in 1:nrow(gseaResult)) print(enrichplot::gseaplot(gseaResult, geneSetID = x, title = paste0(toupper(gseaResult@result$Description[x]), "\nNES = ", sprintf("%.2f", gseaResult@result$NES[x]), " // Adj.p = ", sprintf(fmt = "%.1E", gseaResult@result$p.adjust[x]))))
      dev.off()
      gseaplot.dir <- paste0(gsea.dir, '/gseaplot_Broad/')
      dir.create(path = gseaplot.dir, recursive = TRUE)
      png(paste0(gseaplot.dir, '/GSEA.gseaplot_Broad_%04d.png'), width = 1024, height = 768)
      for (x in 1:nrow(gseaResult)) print(enrichplot::gseaplot2(gseaResult, geneSetID = x, title = paste0(toupper(gseaResult@result$Description[x]), "\nNES = ", sprintf("%.2f", gseaResult@result$NES[x]), " // Adj.p = ", sprintf(fmt = "%.1E", gseaResult@result$p.adjust[x]))))
      dev.off()
    }
    ## KEGG pathview
    if (grepl(pattern = 'gseKEGG', x = gseaResult@setType)) {
      kegg.sp <- tolower(paste0(substr(unlist(strsplit(gseaResult@organism, ' ')), c(1, 1), c(1,2)), collapse = ''))
      kegg.dir <- paste0(gsea.dir, '/KEGG_pathview/')
      dir.create(path = kegg.dir, recursive = TRUE)
      for (x in which(gsea.sig.tf)) {
        kegg.id <- gseaResult@result$ID[x]
        library(pathview)
        pathview::pathview(gene.data = gseaResult@geneList, pathway.id = kegg.id, species = kegg.sp, limit = list(gene=2, cpd=1), low=c("red", "blue"), high = c("cornflowerblue", "yellow"), kegg.dir = kegg.dir)
        file.rename(from = paste0(kegg.id, '.pathview.png'), to = paste0(kegg.dir, '/', kegg.id, '.pathview.png'))
      }
    }
    if (!is.null(cnetplot)) {
      ## Cnetplot (needs readable)
      p <- try(enrichplot::cnetplot(gseaResult_readable, showCategory = cnetplot, circular = FALSE, colorEdge = TRUE), silent = TRUE)
      if(!is(p, class2 = 'try-error')) {
        png(paste0(gsea.dir, '/GSEA.cnetplot_', cnetplot, '.png'), width = 2000, height = 2000)
        print(p + ggplot2::scale_fill_continuous(guide = ggplot2::guide_legend()) + ggplot2::theme(legend.position="bottom"))
        dev.off()
      } else message(paste0('Error captured for cnetplot : \n"', p, '"'))
    }
    if (!is.null(emapplot) & length(which(gsea.sig.tf)) > 1) {
      ## Emapplot (needs readable)
      p <- try(enrichplot::emapplot(enrichplot::pairwise_termsim(gseaResult_readable, showCategory = emapplot), showCategory = emapplot, layout = 'kk'), silent = TRUE)
      if(!is(p, class2 = 'try-error')) {
        png(paste0(gsea.dir, '/GSEA.emapplot_', emapplot, '.png'), width = 1500, height = 1500)
        print(p + ggplot2::scale_fill_continuous(guide = ggplot2::guide_legend()) + ggplot2::theme(legend.position="bottom"))
        dev.off()
      } else message(paste0('Error captured for emapplot : \n"', p, '"'))
    }
  } else message("No enriched term found.")
}

### FUNCTION TO PERFORM ORA
### . 'gene' : character vector : a vector of EntrezIDs corresponding to the topN signature to assess for ORA
### . 'species' : character ; a species name, as in the 'species_name' column of msigdbr::msigdbr_species()
### . 'category' : character ; an MSigDB category, as in the 'gs_cat' column of msigdbr::msigdbr_collections()
### . 'subcategory' : character ; an MSigDB category, as in the 'gs_subcat' column of msigdbr::msigdbr_collections()
### . 'gene2Symbol' : character vector ; a named vector of EntrezIDs, with corresponding Symbol or EnsemblID as names
### . 'pvalueCutoff' : numeric ; minimum FDR-adjusted p-value for signficantly enriched terms (see clusterProfiler::GSEA)
### . 'minGSSize' : integer ; minimal size of each geneSet for analyzing (see clusterProfiler::GSEA)
### . '...' : any other parameter to see clusterProfiler::GSEA
ora.run <- function(gene = NULL, func.name = 'clusterProfiler::enricher', species = 'Homo sapiens', t2g = NULL, t2g.name = NULL, gene2Symbol = NULL, pvalueCutoff = 5E-02, minGSSize = 10, ...) {

  ## Checks
  if (any(!is.character(gene))) stop("'gene' should be a character vector.")
  if (!is.character(func.name)) stop("'func.name' should be a character.")
  if (length(gene) < minGSSize) stop("Length of 'gene' should be at least the value of 'minGSSize'.")
  if (!is.character(species)) stop("'species' should be a character.")
  valid.species <- msigdbr::msigdbr_species()$species_name
  if (!(species %in% valid.species)) stop(paste0("Unsupported 'species'. Expecting one of : '", paste(valid.species, collapse = "', '"), '.'))
  if (!is.null(t2g) & is.null(t2g.name)) stop("'t2g' provided but no 't2g.name'.")
  if (!is.null(t2g.name) & is.null(t2g)) stop("'t2g.name' provided but no 't2g'.")
  if (!is.null(gene2Symbol)) {
    if (any(!is.character(gene2Symbol))) stop("'gene2Symbol' shoud be a character vector, named with EntrezIDs.")
    if (length(names(gene2Symbol)) == 0) stop("'gene2Symbol' shoud be a character vector, named with EntrezIDs.")
  }
  if (!is.numeric(pvalueCutoff)) stop("'pvalueCutoff' should be a positive, < 1, numeric")
  if (!(pvalueCutoff >= 0)) stop("'pvalueCutoff' should be a positive, < 1, numeric")
  if (!(pvalueCutoff <= 1)) stop("'pvalueCutoff' should be a positive, < 1, numeric")
  if (!is.numeric(minGSSize)) stop("'minGSSize' should be an integer.")
  if (!(minGSSize > 0)) stop("'minGSSize' should a non-null positive integer.")
  if (minGSSize < 10) warning("'minGSSize' < 10 : expect no result !")
  
  ## Loading the requested function
  func.split <- unlist(strsplit(func.name, '::'))
  or.function <- base::get(func.split[2], envir = loadNamespace(func.split[1]))

  if (func.name == 'clusterProfiler::enricher') {
    if (is.null(t2g)) stop("When calling this function with func.name = 'clusterProfiler::enricher', a value is required for 't2g'")
    if (is.null(t2g.name)) stop("When calling this function with func.name = 'clusterProfiler::enricher', a value is required for 't2g.name'")
    ## Functions requiring a TERM2GENE (msigdbr, CellMarkers, ...)
    ora.res <- clusterProfiler::enricher(gene = gene, TERM2GENE = t2g, pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, universe = AnnotationDbi::mappedkeys(eval(parse(text = paste0(msigdbr2org(species), 'ACCNUM')))), ...)
    ora.res@ontology <- paste(c(func.name, t2g.name), collapse = '_')
  } else if (func.name == 'meshes::gseMeSH') {
    ## MeSH (requires additional 'MeSHDb', 'database' and 'category' parameters)
    mesh.sp <- paste0(c('MeSH.', substr(unlist(strsplit(species, ' ')), c(1, 1), c(1,2)), '.eg.db'), collapse = '')
    ora.res <- clusterProfiler::enricher(gene = gene, MeSHDb = mesh.sp, pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, universe = AnnotationDbi::mappedkeys(eval(parse(text = paste0(msigdbr2org(species), 'ACCNUM')))), ...)
    ora.res@ontology <- func.name
  } else {
    if ('kegg' %in% tolower(func.name)) {
      ## KEGG / KEGGM (requires a custom species name in 'organism' parameter)
      kegg.sp <- tolower(paste0(substr(unlist(strsplit(species, ' ')), c(1, 1), c(1,2)), collapse = ''))
      ora.res <- or.function(gene = gene, organism = kegg.sp, pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, universe = AnnotationDbi::mappedkeys(eval(parse(text = paste0(msigdbr2org(species), 'ACCNUM')))), ...)
    } else if (species != 'Homo sapiens') stop(paste0("Function '", func.name, "' can only be used with the 'Homo sapiens' species.")) else {
      ## OTHER (generic functions without TERM2GENE or custom parameters, like those in DOSE package)
      ora.res <- or.function(gene = gene, pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, universe = AnnotationDbi::mappedkeys(eval(parse(text = paste0(msigdbr2org(species), 'ACCNUM')))), ...)
    }
    ora.res@ontology <- func.name
  }

  ## Adding some useful metadata
  ora.res@organism <- species
  if(!is.null(gene2Symbol)) {
    ora.res@gene2Symbol <- gene2Symbol
    ora.res@keytype <- 'ENTREZID'
  }
  return(ora.res)
}

### FUNCTION TO OUTPUT TABLE AND PLOT RESULTS FROM ora.run()
### . 'enrichResults' : enrichResults object ; an output from the ora.run() function
### . 'comp.name' : character ; The name of the differential expression comparison which resulted in the 'gene' input to ora.run(). Actually, only used in plots' title.
### . 'out.dir' : character ; directory to output plots and tables
### . 'heatplot' : integer of NULL ; if integer perform an enrichplot::dotplot() limited to a number of term corresponding to the given integer value. If null, no plot is performed.
### . 'geneList' : numeric vector : a named vector of decreasing values, with EntrezIDs as names (the one object that is used as the 'geneList' input for the gsea.run() function). Required if 'heatplot' is non-NULL, and for KEGG pathview.
### . 'dotplot' : integer of NULL ; if integer perform an enrichplot::dotplot() limited to a number of term corresponding to the given integer value. If null, no plot is performed.
### . 'barplot' : integer of NULL ; if integer perform an enrichplot::barplot() limited to a number of term corresponding to the given integer value. If null, no plot is performed.
### . 'cnetplot' : integer of NULL ; if integer perform an enrichplot::cnetplot() limited to a number of term corresponding to the given integer value. If null, no plot is performed.
### . 'emapplot' : integer of NULL ; if integer perform an enrichplot::emapplot() limited to a number of term corresponding to the given integer value. If null, no plot is performed.
ora.output <- function(enrichResult = NULL, comp.name = 'TEST', out.dir = getwd(), heatplot = 100, geneList = NULL, dotplot = 100, barplot = 100, cnetplot = 10, emapplot = 10) {

  ## Significant terms ?
  ora.sig.tf <- enrichResult@result$p.adjust < enrichResult@pvalueCutoff

  if (any(ora.sig.tf)) {

    ## Creating ORA output dir
    ora.dir <- paste(c(out.dir, paste0('ORA_adjp.', enrichResult@pvalueCutoff), enrichResult@ontology), collapse = '/')
    dir.create(path = ora.dir, recursive = TRUE)
    
    ## Dumping ORA results
    saveRDS(enrichResult, file = paste0(ora.dir, '/ORA.results.RDS'), compress = 'bzip2')
    
    ## Converting to readable symbols
    enrichResult_readable <- DOSE::setReadable(enrichResult, OrgDb = paste0(msigdbr2org(gseaResult@organism), '.db'))
    
    ## Writing readable table
    write.table(enrichResult_readable@result, file = paste0(ora.dir, '/ORA.results_readable.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
    
    ## Heatplot (needs readable)
    if (!is.null(heatplot) & !is.null(geneList)) {
      hp <- enrichplot::heatplot(enrichResult_readable, showCategory = heatplot, foldChange = geneList)
      png(paste0(ora.dir, '/ORA.heatplot.png'), width = (5 * nrow(hp$data)) + 500, height = (min(length(which(ora.sig.tf)), heatplot) * 12) + 150)
      print(hp)
      dev.off()
    }
    
    ## Dotplot // Barplot
    complot <- c()
    if (!is.null(dotplot)) complot <- c(complot, setNames(dotplot, 'clusterProfiler::dotplot'))
    if (!is.null(barplot)) complot <- c(complot, setNames(barplot, 'enrichplot:::barplot.enrichResult'))
    for (cx in seq_along(complot)) {
      func.split <- unlist(strsplit(names(complot)[cx], ':+'))
      plot.function <- base::get(func.split[2], envir = loadNamespace(func.split[1]))
      enrichResult_bak <- enrichResult
      enrichResult@result$Description <- gsub(pattern = '_', ' ', enrichResult@result$Description, fixed = TRUE)
      png(paste0(ora.dir, '/ORA.', func.split[2], '.png'), width = 1000, height = (min(length(which(ora.sig.tf)), unname(complot[cx])) * 20) + 300)
      print(plot.function(enrichResult, showCategory = unname(complot[cx]), title = paste(c(comp.name, enrichResult@ontology, 'ORA'), collapse = "\n")) + ggplot2::scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=100)))
      dev.off()
      enrichResult <- enrichResult_bak
      rm(enrichResult_bak)
    }
    ## KEGG pathview
    if (grepl(pattern = 'enrichKEGG', x = enrichResult@ontology) & !is.null(geneList)) {
      kegg.sp <- tolower(paste0(substr(unlist(strsplit(enrichResult@organism, ' ')), c(1, 1), c(1,2)), collapse = ''))
      kegg.dir <- paste0(ora.dir, '/KEGG_pathview/')
      dir.create(path = kegg.dir, recursive = TRUE)
      for (x in which(enrichResult@result$p.adjust < enrichResult@pvalueCutoff)) {
        kegg.id <- enrichResult@result$ID[x]
        ## Custom geneList limited to our topN signature
        cur.geneList <- geneList[enrichResult@gene2Symbol]
        library(pathview)
        suppressMessages(pathview::pathview(gene.data = cur.geneList, pathway.id = kegg.id, species = kegg.sp, limit = list(gene=2, cpd=1), low=c("red", "blue"), high = c("cornflowerblue", "yellow"), kegg.dir = kegg.dir))
        file.rename(from = paste0(kegg.id, '.pathview.png'), to = paste0(kegg.dir, '/', kegg.id, '.pathview.png'))
      }
    }
    if (!is.null(cnetplot)) {
      ## Cnetplot (need readable)
      p <- try(enrichplot::cnetplot(enrichResult_readable, showCategory = cnetplot, circular = FALSE, colorEdge = TRUE), silent = TRUE)
      if(!is(p, class2 = 'try-error')) {
        png(paste0(ora.dir, '/GSEA.cnetplot_', cnetplot, '.png'), width = 2000, height = 2000)
        print(p + ggplot2::scale_fill_continuous(guide = ggplot2::guide_legend()) + ggplot2::theme(legend.position="bottom"))
        dev.off()
      } else message(paste0('Error captured for cnetplot : \n"', p, '"'))
    }
    if (!is.null(emapplot)  & length(which(ora.sig.tf)) > 1) {
      ## Emapplot (need readable)
      p <- try(enrichplot::emapplot(enrichplot::pairwise_termsim(enrichResult_readable, showCategory = emapplot), showCategory = emapplot, layout = 'kk'), silent = TRUE)
      if(!is(p, class2 = 'try-error')) {
        png(paste0(ora.dir, '/GSEA.emapplot_', emapplot, '.png'), width = 2000, height = 2000)
        print(p + ggplot2::scale_fill_continuous(guide = ggplot2::guide_legend()) + ggplot2::theme(legend.position="bottom"))
        dev.off()
      } else message(paste0('Error captured for emapplot : \n"', p, '"'))
      
    }
  } else message("No enriched term found.")
}



##################
#### EXAMPLES ####
##################

## Setting variables
defile <- '/home/job/WORKSPACE/B21002_DELE_01/B21002_DELE_01_dge/DGE/GSEAapp/Macrophage_type_compairing_Proinf_vs_Basal/Macrophage_type_compairing_Proinf_vs_Basal_complete.tsv'
out.dir <- dirname(defile)
comp.name <- basename(dirname(defile))
species <- 'Homo sapiens'
my.seed <- 1337
de.min.p <- 5E-02
enr.min.p <- 5E-02
enr.min.genes <- 10
topN <- 100

## PREPARING INPUT (from a '*_complete.tsv' output table from our 'rna-salmon-deseq2' pipeline)
enr.inputs <- pipe2enr(deseq2.res.file = defile, species = species, geneid.colname = 'Gene_Name', geneid.type = 'SYMBOL', value.colname = 'stat_change', topN.max = topN, topN.order.colname = 'Adjusted_PValue', topN.order.decreasing = FALSE, topN.cutoff = de.min.p, topN.keep.operator = '<', header = TRUE, sep = "\t", as.is = TRUE)

## MSIGDB
### Query any (and in this example, all) MsigDB banks available in the 'msigdbr' package. See http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
### Get available species / collections
msigdb.collections <- as.data.frame(msigdbr::msigdbr_collections())
msigdb.species <- as.data.frame(msigdbr::msigdbr_species())

### GSEA
for (x in seq_len(nrow(msigdb.collections))) {
  my.collec <- unlist(msigdb.collections[x, c("gs_cat", "gs_subcat"), drop = TRUE])
  ## Import the TERM2GENE object corresponding to the desired category/subcategory combo
  my.t2g <- msigdb_to_t2g(species = species, category = my.collec[1], subcategory = my.collec[2])
  my.t2g.name <- unname(if(my.collec[2] == '') my.collec[1] else paste(my.collec, collapse = "_"))
  ## Run the GSEA
  my.gsea.res <- gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = 'clusterProfiler::GSEA', t2g = my.t2g, t2g.name = my.t2g.name, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, nPerm = nPerm, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes)
  ## Generate plots / outputs
  gsea.output(gseaResult = my.gsea.res, out.dir = out.dir, comp.name = comp.name)
}

### ORA
for (x in seq_len(nrow(msigdb.collections))) {
  my.collec <- unlist(msigdb.collections[x, c("gs_cat", "gs_subcat"), drop = TRUE])
  ## Import the TERM2GENE object corresponding to the desired category/subcategory combo
  my.t2g <- msigdb_to_t2g(species = species, category = my.collec[1], subcategory = my.collec[2])
  my.t2g.name <- ifelse(my.collec[2] == '', my.collec[1], paste(my.collec, collapse = "_"))
  ## Run the ORA
  my.ora.res <- ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = 'clusterProfiler::enricher', t2g = my.t2g, t2g.name = my.t2g.name, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes)
  ## Generate plots / outputs
  ora.output(enrichResult = my.ora.res, out.dir = out.dir, comp.name = comp.name, geneList = enr.inputs$gsea.genevec)
}

## DO (Disease Ontology), NCG (Network of Cancer Genes), DGN (DisGeNET)
### WARNING : ONLY FOR 'Homo sapiens' !

### GSEA
for (x in c('DOSE::gseDO', 'DOSE::gseNCG', 'DOSE::gseDGN')) {
  my.gsea.res <- gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes)
  gsea.output(gseaResult = my.gsea.res, out.dir = out.dir, comp.name = comp.name)
}

### ORA
for (x in c('DOSE::enrichDO', 'DOSE::enrichNCG', 'DOSE::enrichDGN')) {
  my.ora.res <- ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes)
  ora.output(enrichResult = my.ora.res, out.dir = out.dir, comp.name = comp.name, geneList = enr.inputs$gsea.genevec)
}

## KEGG/MKEGG
### NOTE1 : It's the same way to call the 'gsea.run' / 'ora.run' as it is for 'DO', 'NCG' or 'DGN', but here it's compatible with many more species than homo sapiens.
### NOTE2 : for this case, additional KEGG pathway plots will be generated.
### NOTE3 : for this case, an internet connexion is required to query the KEGG website.

### GSEA
for (x in c('clusterProfiler::gseKEGG', 'clusterProfiler::gseMKEGG')) {
  my.gsea.res <- gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes)
  gsea.output(gseaResult = my.gsea.res, out.dir = out.dir, comp.name = comp.name)
}

### ORA
for (x in c('clusterProfiler::enrichKEGG', 'clusterProfiler::enrichMKEGG')) {
  my.ora.res <- ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes)
  ora.output(enrichResult = my.ora.res, out.dir = out.dir, comp.name = comp.name, geneList = enr.inputs$gsea.genevec)
}

## MESH (WARNING : VERY VERY SLOW as big DBs, 3 sources, 16 categories !)
### Requires additional parameters :
### . 'MeSHDb' : character ; name of a MeSH [NO : AUTO FROM SPECIES NAME]
### . 'database' : character ; MeSH source type (can be 'gendoo' = text-mining, 'gene2pubmed' = manual curation by NCBI team, 'RBBH' = sequence homology with BLASTP search @ E-value < 1E-50)
### . 'category' : character ; name of a MeSH category sub-db.
### NOTE : see https://yulab-smu.top/biomedical-knowledge-mining-book/meshes-semantic-similarity.html
mesh.func.name <- 'meshes::gseMeSH'
mesh.dbs <- c('gendoo', 'gene2pubmed', 'RBBH')
mesh.categories <- toupper(letters[-c(15:21,23:25)])

### GSEA
for (y in mesh.dbs) {
  for (x in mesh.categories) {
    my.gsea.res <- gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = mesh.func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes, database = y, category = x)
    ## Little hack specific to MeSH results (as I was not able to get the value of extra parameters 'database' and 'category' from within the 'gsea.run()' function)
    my.gsea.res@setType <- paste(c(mesh.func.name, y, x), collapse = '_')
    gsea.output(gseaResult = my.gsea.res, out.dir = out.dir, comp.name = comp.name)
  }
}

### ORA
for (i in mesh.dbs) {
  for (x in mesh.categories) {
    my.ora.res <- ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = mesh.func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes, database = y, category = x)
    ## Little hack specific to MeSH results (as I was not able to get the value of extra parameters 'database' and 'category' from within the 'gsea.run()' function)
    my.ora.res@ontology <- paste(c(mesh.func.name, database, category), collapse = '_')
    ora.output(enrichResult = my.ora.res, out.dir = out.dir, comp.name = comp.name, geneList = enr.inputs$gsea.genevec)
  }
}

## CELLMARKERS
### Assess cell types from an online table
### NOTE : It's the same way to call the 'gsea.run' / 'ora.run' functions as for MSIGDB, but with a single bank (so, no loop).

### Import the CellMarkers bank
library(tidyverse)
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>% tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% dplyr::select(cellMarker, geneID) %>% dplyr::mutate(geneID = strsplit(geneID, ', '))

### GSEA
my.gsea.res <- gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = 'clusterProfiler::GSEA', t2g = cell_markers, t2g.name = 'CellMarkers', gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes)
gsea.output(gseaResult = my.gsea.res, out.dir = out.dir, comp.name = comp.name)

#### ORA
my.ora.res <- ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = 'clusterProfiler::enricher', t2g = cell_markers, t2g.name = 'CellMarkers', gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes)
ora.output(enrichResult = my.ora.res, out.dir = out.dir, comp.name = comp.name, geneList = enr.inputs$gsea.genevec)




##################
#### TESTZONE ####
##################

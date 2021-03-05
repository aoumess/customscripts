################
#### HEADER ####
################

## A Wrapper to automate GSEA/ORA functional annotation and analysis, based on clusterProfiler
## PACKAGES :
### REQUIRED :
#### CRAN : 'tidyverse'
#### Bioconductor : 'clusterProfiler' (will import other required packages as dependencies, like 'DOSE' and 'enrichplot'), 'org.Xx.eg.db' (species-specific, like 'org.Hs.eg.db')
### SUGGESTED :
#### Bioconductor : 'msigdbr' (to assess MSigDB databases), 'meshes' (to assess MeSH databases), 'MeSH.Xxx.eg.db' (species-specific, like 'MeSH.Hsa.eg.db', required to query MeSH databases), 'pathview' (to plot KEGG pathways if 'clusterProfiler::enrichKEGG' or 'clusterProfiler::gseKEGG' functions are used)




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

### FUNCTION TO PERFORM GSEA
### . 'geneList' : numeric vector : a named vector of decreasing values, with EntrezIDs as names
### . 'species' : character ; a species name, as in the 'species_name' column of msigdbr::msigdbr_species()
### . 'category' : character ; an MSigDB category, as in the 'gs_cat' column of msigdbr::msigdbr_collections()
### . 'subcategory' : character ; an MSigDB category, as in the 'gs_subcat' column of msigdbr::msigdbr_collections()
### . 'gene2Symbol' : character vector ; a named vector of EntrezIDs, with corresponding Symbol or EnsemblID as names
### . 'seed' : integer ; the random number generation seed
### . 'nPerm' : integer ; number of permutations for GSEA (see clusterProfiler::GSEA)
### . 'pvalueCutoff' : numeric ; minimum FDR-adjusted p-value for signficantly enriched terms (see clusterProfiler::GSEA)
### . 'minGSSize' : integer ; minimal size of each geneSet for analyzing (see clusterProfiler::GSEA)
### . '...' : any other parameter to pass to 'func.name'
### NOTE : For input ('geneList'), ALWAYS use all available genes (ie, not limited to significant differentialy expressed genes), as intended for GSEA analysis. NEVER use a selection of genes (results will be corrupt). Please also be careful that you may obtain significant GSEA results from a list of values corresponding to NO significant genes (as GSEA only needs the order of these genes). In this case you may talk about "tendencies" of enrichment in your differential expression results.
gsea.run <- function(geneList = NULL, func.name = 'clusterProfiler::GSEA', species = 'Homo sapiens', t2g = NULL, t2g.name = NULL, gene2Symbol = NULL, seed = 1337, nPerm = 1E+05, pvalueCutoff = 5E-02, minGSSize = 10, ...) {

  ## Loading the requested function
  func.split <- unlist(strsplit(func.name, '::'))
  gse.function <- base::get(func.split[2], envir = loadNamespace(func.split[1]))

  if (func.name == 'clusterProfiler::GSEA') {
    if (is.null(t2g)) stop("When calling this function with func.name = 'clusterProfiler::GSEA', a value is required for 't2g'")
    if (is.null(t2g.name)) stop("When calling this function with func.name = 'clusterProfiler::enricher', a value is required for 't2g.name'")
    
    ## Functions requiring a TERM2GENE (msigdbr, CellMarkers, ...)
    gsea.res <- gse.function(geneList = geneList, TERM2GENE = t2g, pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, seed = seed, nPerm = nPerm, ...)
    gsea.res@setType <- paste(c(func.name, t2g.name), collapse = '_')
  } else if (func.name == 'meshes::gseMeSH') {
    ## MeSH (requires additional 'MeSHDb', 'database' and 'category' parameters)
    mesh.sp <- paste0(c('MeSH.', substr(unlist(strsplit(species, ' ')), c(1, 1), c(1,2)), '.eg.db'), collapse = '')
    gsea.res <- gse.function(geneList = geneList, MeSHDb = mesh.sp, pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, seed = seed, nPerm = nPerm, ...)
    gsea.res@setType <- paste(c(func.name, database, category), collapse = '_')
  } else {
    if ('kegg' %in% tolower(func.name)) {
      ## KEGG / KEGGM (requires a custom species name in 'organism' parameter)
      kegg.sp <- tolower(paste0(substr(unlist(strsplit(species, ' ')), c(1, 1), c(1,2)), collapse = ''))
      gsea.res <- gse.function(geneList = geneList, organism = kegg.sp, pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, seed = seed, nPerm = nPerm, ...)
    } else if (species != 'Homo sapiens') stop(paste0("Function '", func.name, "' can only be used with the 'Homo sapiens' species.")) else {
      ## OTHER (generic functions without TERM2GENE or custom parameters, like those in DOSE package)
      gsea.res <- gse.function(geneList = geneList, pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, seed = seed, nPerm = nPerm, ...)
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
### . 'heatplot.do' : logic ; perform an enrichplot::heatplot()
### . 'dotplot.do' : logic ; perform an enrichplot::dotplot()
### . 'plot.maxterms' : integer ; Maximum number of displayed enriched terms
### . 'out.dir' : character ; directory to output plots and tables
gsea.output <- function(gseaResult = NULL, comp.name = 'TEST', heatplot.do = TRUE, dotplot.do = TRUE, plot.maxterms = 100, barplot = TRUE, gsea.plot = TRUE, ridgeplot.do = TRUE, cnetplot = 10, emapplot = 10, out.dir = getwd()) {

  ## Significant terms ?
  gsea.sig.tf <- gseaResult@result$p.adjust < gseaResult@params$pvalueCutoff

  if (any(gsea.sig.tf)) {

    ## Creating GSEA output dir
    gsea.dir <- paste(c(out.dir, paste0('GSEA_adjp.', gseaResult@params$pvalueCutoff, '_nPerm.', gseaResult@params$nPerm), gseaResult@setType), collapse = '/')
    dir.create(path = gsea.dir, recursive = TRUE)

    ## Dumping GSEA results
    saveRDS(gseaResult, file = paste0(gsea.dir, '/GSEA.results.RDS'), compress = 'bzip2')
    
    ## Converting to readable symbols
    gseaResult_readable <- DOSE::setReadable(gseaResult, OrgDb = paste0(msigdbr2org(gseaResult@organism), '.db'))
    ## Writing readable table
    write.table(gseaResult_readable@result, file = paste0(gsea.dir, '/GSEA.results_readable.txt'), quote = FALSE, sep = "\t", row.names = FALSE)

    ## Heatplot (needs readable)
    if (heatplot.do) {
      hp <- clusterProfiler::heatplot(gseaResult_readable, foldChange = gseaResult@geneList)
      png(paste0(gsea.dir, '/GSEA.heatplot.png'), width = (5 * nrow(hp$data)) + 500, height = (min(length(which(gsea.sig.tf)), plot.maxterms) * 12) + 150)
      print(hp)
      dev.off()
    }
    
    ## Dotplot // Barplot
    complot <- c()
    if (dotplot.do) complot <- c(complot, 'clusterProfiler::dotplot')
    if (barplot.do) complot <- c(complot, 'enrichplot:::barplot.enrichResult')
    for (cx in complot) {
      func.split <- unlist(strsplit(cx, ':+'))
      plot.function <- base::get(func.split[2], envir = loadNamespace(func.split[1]))
      gseaBak <- gseaResult
      gseaResult@result$Description <- gsub(pattern = '_', ' ', gseaResult@result$Description, fixed = TRUE)
      png(paste0(gsea.dir, '/GSEA.', func.split[2], '.png'), width = 1000, height = (min(length(which(gsea.sig.tf)), plot.maxterms) * 20) + 300)
      print(plot.function(gseaResult, showCategory = plot.maxterms, title = paste(c(comp.name, gseaResult@setType, 'GSEA'), collapse = "\n"), split = '.sign') + ggplot2::facet_grid(~ .sign) + ggplot2::scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=100)))
      dev.off()
      gseaResult <- gseaBak
      rm(gseaBak)
    }
    
    if (ridgeplot.do) {
      ## Ridgeplot
      gseaBak <- gseaResult
      gseaResult@result$Description <- gsub(pattern = '_', ' ', gseaResult@result$Description, fixed = TRUE)
      rip <- enrichplot::ridgeplot(gseaResult, showCategory = plot.maxterms)
      png(paste0(gsea.dir, '/GSEA.ridgeplot.png'), width = 1000, height = (min(length(which(gsea.sig.tf)), plot.maxterms) * 20) + 300)
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
      png(paste0(gsea.dir, '/GSEA.cnetplot_', cnetplot, '.png'), width = 2000, height = 2000)
      print(enrichplot::cnetplot(gseaResult_readable, showCategory = cnetplot, circular = FALSE, colorEdge = TRUE) + ggplot2::scale_fill_continuous(guide = guide_legend()) + ggplot2::theme(legend.position="bottom"))
      dev.off()
    }
    if (!is.null(emapplot) & length(which(gsea.sig.tf)) > 1) {
      ## Emapplot (needs readable)
      png(paste0(gsea.dir, '/GSEA.emapplot_', emapplot, '.png'), width = 1500, height = 1500)
      print(enrichplot::emapplot(gseaResult_readable, showCategory = emapplot, layout = 'kk') + ggplot2::scale_fill_continuous(guide = guide_legend()) + ggplot2::theme(legend.position="bottom"))
      dev.off()
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
### . 'setReadable' : logic ; output results with Symbol (or EnsemblIDs, depending on what was provided to the 'gene2Symbol' parameter) instead of EntrezIDs
### . '...' : any other parameter to see clusterProfiler::GSEA
ora.run <- function(gene = NULL, func.name = 'clusterProfiler::enricher', species = 'Homo sapiens', t2g = NULL, t2g.name = NULL, gene2Symbol = NULL, pvalueCutoff = 5E-02, minGSSize = 10, setReadable = TRUE, ...) {

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
    ora.res@ontology <- paste(c(func.name, database, category), collapse = '_')
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
  if (setReadable) ora.res <- DOSE::setReadable(ora.res, OrgDb = paste0(msigdbr2org(species), '.db'))
  return(ora.res)
}

### FUNCTION TO OUTPUT TABLE AND PLOT RESULTS FROM ora.run()
### . 'enrichResults' : enrichResults object ; an output from the ora.run() function
### . 'comp.name' : character ; The name of the differential expression comparison which resulted in the 'gene' input to ora.run(). Actually, only used in plots' title.
### . 'heatplot.do' : logic ; perform an enrichplot::heatplot(). If TRUE, a 'geneList' is required.
### . 'geneList' : numeric vector : a named vector of decreasing values, with EntrezIDs as names (the one object that is used as the 'geneList' input for the gsea.run() function). Required if 'heatplot.do' set to TRUE and for KEGG pathview.
### . 'dotplot.do' : logic ; perform an enrichplot::dotplot()
### . 'plot.maxterms' : integer ; Maximum number of displayed enriched terms
### . 'out.dir' : character ; directory to output plots and tables
ora.output <- function(enrichResult = NULL, comp.name = 'TEST', heatplot.do = TRUE, geneList = NULL, dotplot.do = TRUE, plot.maxterms = 100, barplot.do = TRUE, cnetplot = 10, emapplot = 10, out.dir = getwd()) {

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
    if (heatplot.do & !is.null(geneList)) {
      hp <- clusterProfiler::heatplot(enrichResult_readable, foldChange = geneList)
      png(paste0(ora.dir, '/ORA.heatplot.png'), width = (5 * nrow(hp$data)) + 500, height = (min(length(which(ora.sig.tf)), plot.maxterms) * 12) + 150)
      print(hp)
      dev.off()
    }
    
    ## Dotplot // Barplot
    complot <- c()
    if (dotplot.do) complot <- c(complot, 'clusterProfiler::dotplot')
    if (barplot.do) complot <- c(complot, 'graphics::barplot')
    for (x in complot) {
      func.split <- unlist(strsplit(x, '::'))
      plot.function <- base::get(func.split[2], envir = loadNamespace(func.split[1]))
      enrichResult_bak <- enrichResult
      enrichResult@result$Description <- gsub(pattern = '_', ' ', enrichResult@result$Description, fixed = TRUE)
      png(paste0(ora.dir, '/ORA.', func.split[2], '.png'), width = 1000, height = (min(length(which(ora.sig.tf)), plot.maxterms) * 20) + 300)
      print(plot.function(enrichResult, showCategory = plot.maxterms, title = paste(c(comp.name, enrichResult@ontology, 'ORA'), collapse = "\n")) + ggplot2::scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=100)))
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
        cur.geneList <- geneList[names(enrichResult@gene2Symbol)]
        library(pathview)
        suppressMessages(pathview::pathview(gene.data = cur.geneList, pathway.id = kegg.id, species = kegg.sp, limit = list(gene=2, cpd=1), low=c("red", "blue"), high = c("cornflowerblue", "yellow"), kegg.dir = kegg.dir))
        file.rename(from = paste0(kegg.id, '.pathview.png'), to = paste0(kegg.dir, '/', kegg.id, '.pathview.png'))
      }
    }
    if (!is.null(cnetplot)) {
      ## Cnetplot (need readable)
      png(paste0(ora.dir, '/GSEA.cnetplot_', cnetplot, '.png'), width = 2000, height = 2000)
      print(enrichplot::cnetplot(enrichResult_readable, showCategory = cnetplot, circular = FALSE, colorEdge = TRUE) + ggplot2::scale_fill_continuous(guide = guide_legend()) + ggplot2::theme(legend.position="bottom"))
      dev.off()
    }
    if (!is.null(emapplot)) {
      ## Cnetplot (need readable)
      png(paste0(ora.dir, '/GSEA.emapplot_', emapplot, '.png'), width = 2000, height = 2000)
      print(enrichplot::emapplot(enrichResult_readable, showCategory = emapplot, layout = 'kk') + ggplot2::scale_fill_continuous(guide = guide_legend()) + ggplot2::theme(legend.position="bottom"))
      dev.off()
    }
  } else message("No enriched term found.")
}




#####################
#### PREPARATION ####
#####################

## Setting variables
species <- 'Homo sapiens'
gsea.do <- TRUE
ora.do <- TRUE
my.seed <- 1337
nPerm <- 1E+05
enr.min.p <- 5E-02
enr.min.genes <- 10
topN <- 100

## Preparing input
source.genetype <- 'SYMBOL' ## ('SYMBOL' or 'ENSEMBLID')
testfile <- '/home/job/WORKSPACE/B21002_DELE_01/B21002_DELE_01_dge/DGE/GSEAapp/Macrophage_type_compairing_Proinf_vs_Basal/Macrophage_type_compairing_Proinf_vs_Basal_complete.tsv'
value.colname <- 'stat_change' ## Here : logFC
testdata <- read.table(testfile, header = TRUE, sep = "\t", as.is = TRUE)
gconv <- as.list(clusterProfiler::bitr(testdata$Gene_Name, fromType=source.genetype, toType='ENTREZID', OrgDb=msigdbr2org(species)))
names(gconv$ENTREZID) <- gconv[[source.genetype]]
names(gconv[[source.genetype]]) <- gconv$ENTREZID
testdata$ENTREZID <- gconv$ENTREZID[as.character(testdata$Gene_Name)]
### Removing entries without EntrezID
testdata <- testdata[!is.na(testdata$ENTREZID),]
### Removing entries with duplicated EntrezID
testdata <- testdata[!testdata$ENTREZID %in% unique(testdata$ENTREZID[duplicated(testdata$ENTREZID)]),]
### GSEA-formatted input
gsea.genevec <- sort(setNames(testdata$stat_change, testdata$ENTREZID), decreasing = TRUE)
### ORA-formatted input
ora.genevec <- testdata[order(testdata$Adjusted_PValue), "ENTREZID"][1:topN]

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
### . '...' : any parameter to read.table() to handle the input file
pipe2gsea <- function(deseq2.res.file = NULL, species = "Homo sapiens", geneid.colname = 'Gene_Name', geneid.type = 'SYMBOL', value.colname = 'stat_change', topN.max = 100, topN.order.colname = 'Adjusted_PValue', topN.order.decreasing = FALSE, topN.keep.limit = 5E-02, topN.keep.operator = '<', ...) {
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
  deres2 <- deres2[eval(parse(text = paste0('deres2[["', topN.order.colname, '"]] ', topN.keep.operator, ' topN.cutoff'))),]
  
  return(list(gsea.genevec = gsea.genevec,
              ora.genevec = setNames(deres[order(deres[[order.by.colname]], decreasing = order.by.decreasing), "ENTREZID"][1:topN], deres[order(deres[[order.by.colname]], decreasing = order.by.decreasing), geneid.colname][1:topN])))
}


##################
#### EXAMPLES ####
##################

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
  my.t2g.name <- ifelse(my.collec[2] == '', my.collec[1], paste(my.collec, collapse = "_"))
  ## Run the GSEA
  my.gsea.res <- gsea.run(geneList = gsea.genevec, species = species, func.name = 'clusterProfiler::GSEA', t2g = my.t2g, t2g.name = my.t2g.name, gene2Symbol = gconv$ENTREZID, seed = my.seed, nPerm = nPerm, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes)
  ## Generate plots / outputs
  gsea.output(gseaResult = my.gsea.res, out.dir = dirname(testfile), comp.name = basename(dirname(testfile)))
}

### ORA
for (x in seq_len(nrow(msigdb.collections))) {
  my.collec <- unlist(msigdb.collections[x, c("gs_cat", "gs_subcat"), drop = TRUE])
  ## Import the TERM2GENE object corresponding to the desired category/subcategory combo
  my.t2g <- msigdb_to_t2g(species = species, category = my.collec[1], subcategory = my.collec[2])
  my.t2g.name <- ifelse(my.collec[2] == '', my.collec[1], paste(my.collec, collapse = "_"))
  ## Run the ORA
  my.ora.res <- ora.run(gene = ora.genevec, species = species, func.name = 'clusterProfiler::enricher', t2g = my.t2g, t2g.name = my.t2g.name, gene2Symbol = gconv$ENTREZID, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes, setReadable = TRUE)
  ## Generate plots / outputs
  ora.output(enrichResult = my.ora.res, out.dir = dirname(testfile), comp.name = basename(dirname(testfile)), geneList = gsea.genevec)
}

## DO (Disease Ontology), NCG (Network of Cancer Genes), DGN (DisGeNET)
### WARNING : ONLY FOR 'Homo sapiens' !

### GSEA
for (x in c('DOSE::gseDO', 'DOSE::gseNCG', 'DOSE::gseDGN')) {
  my.gsea.res <- gsea.run(geneList = gsea.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = gconv$ENTREZID, seed = my.seed, nPerm = nPerm, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes, setReadable = TRUE)
  gsea.output(gseaResult = my.gsea.res, out.dir = dirname(testfile), comp.name = basename(dirname(testfile)))
}

### ORA
for (x in c('DOSE::enrichDO', 'DOSE::enrichNCG', 'DOSE::enrichDGN')) {
  my.ora.res <- ora.run(gene = ora.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = gconv$ENTREZID, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes, setReadable = TRUE)
  ora.output(enrichResult = my.ora.res, out.dir = dirname(testfile), comp.name = basename(dirname(testfile)), geneList = gsea.genevec)
}

## KEGG/MKEGG
### NOTE1 : It's the same way to call the 'gsea.run' / 'ora.run' as it is for 'DO', 'NCG' or 'DGN', but here it's compatible with many more species than homo sapiens.
### NOTE2 : for this case, additional KEGG pathway plots will be generated.
### NOTE3 : for this case, an internet connexion is required to query the KEGG website.

### GSEA
for (x in c('clusterProfiler::gseKEGG', 'clusterProfiler::gseMKEGG')) {
  my.gsea.res <- gsea.run(geneList = gsea.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = gconv$ENTREZID, seed = my.seed, nPerm = nPerm, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes, setReadable = TRUE)
  gsea.output(gseaResult = my.gsea.res, out.dir = dirname(testfile), comp.name = basename(dirname(testfile)))
}

### ORA
for (x in c('clusterProfiler::enrichKEGG', 'clusterProfiler::enrichMKEGG')) {
  my.ora.res <- ora.run(gene = ora.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = gconv$ENTREZID, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes, setReadable = TRUE)
  ora.output(enrichResult = my.ora.res, out.dir = dirname(testfile), comp.name = basename(dirname(testfile)), geneList = gsea.genevec)
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
for (i in mesh.dbs) {
  for (x in mesh.categories) {
    my.gsea.res <- gsea.run(geneList = gsea.genevec, species = species, func.name = mesh.func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = gconv$ENTREZID, seed = my.seed, nPerm = nPerm, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes, setReadable = TRUE, database = y, category = x)
    gsea.output(gseaResult = my.gsea.res, out.dir = dirname(testfile), comp.name = basename(dirname(testfile)))
  }
}

### ORA
for (i in mesh.dbs) {
  for (x in mesh.categories) {
    my.ora.res <- ora.run(gene = ora.genevec, species = species, func.name = mesh.func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = gconv$ENTREZID, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes, setReadable = TRUE, database = y, category = x)
    ora.output(enrichResult = my.ora.res, out.dir = dirname(testfile), comp.name = basename(dirname(testfile)), geneList = gsea.genevec)
  }
}

## CELLMARKERS
### Assess cell types from an online table
### NOTE : It's the same way to call the 'gsea.run' / 'ora.run' functions as for MSIGDB, but with a single bank (so, no loop).

### Import the CellMarkers bank
library(tidyverse)
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>% tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% dplyr::select(cellMarker, geneID) %>% dplyr::mutate(geneID = strsplit(geneID, ', '))

### GSEA
my.gsea.res <- gsea.run(geneList = gsea.genevec, species = species, func.name = 'clusterProfiler::GSEA', t2g = cell_markers, t2g.name = 'CellMarkers', gene2Symbol = gconv$ENTREZID, seed = my.seed, nPerm = nPerm, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes, setReadable = TRUE)
gsea.output(gseaResult = my.gsea.res, out.dir = dirname(testfile), comp.name = basename(dirname(testfile)))

#### ORA
my.ora.res <- ora.run(gene = ora.genevec, species = species, func.name = 'clusterProfiler::enricher', t2g = cell_markers, t2g.name = 'CellMarkers', gene2Symbol = gconv$ENTREZID, pvalueCutoff = enr.min.p, minGSSize = enr.min.genes, setReadable = TRUE)
ora.output(enrichResult = my.ora.res, out.dir = dirname(testfile), comp.name = basename(dirname(testfile)), geneList = gsea.genevec)




##################
#### TESTZONE ####
##################

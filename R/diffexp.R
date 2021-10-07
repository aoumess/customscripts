## Perform DE analysis & functional enrichment for contrasts in pairs
## exp.mat                matrix(integer)     Sample x gene raw count matrix
## annot.df               data.frame          Sample annotations
## samples.colname        character           Name of the annot.df column to identify samples
## batch.colnames         vector(character)   Name(s) of the annot.df columns to consider sources of batch effect (sequencing run, etc) **MAX LENGTH = 2 (as supported by limma::removeBatchEffect)**
## factorname             vector(character)   Name(s) of the annot.df columns to consider factors on which computing the differential analysis
## invert.levels          logical             If TRUE, inverts the default order of levels (if the default one does not fit your biological preferences)
## adjp.max               0<numeric<1         BH FDR-adjusted p-value cut-off to consider differential genes as significant
## lfc.min                numeric+            Minimal logFoldChange to consider differential genes
## outdir                 character           Path of the output directory
## samples.dist.method    character           Name of the distance method to use for the hierarchical clustering of samples
## samples.hclust.method  character           Name of the aggregation method to use for the hierarchical clustering of samples
## genes.dist.method      character           Name of the distance method to use for the hierarchical clustering of genes
## genes.hclust.method    character           Name of the aggregation method to use for the hierarchical clustering of genes
## msigdb.do              logical, logical    Use the MSigDb dataset collection to perform GSEA/ORA
## go.do                  logical, logical    Use the GeneOntology dataset to perform GSEA/ORA
## do.do                  logical, logical    Use the DiseaseOntology + CancerGeneNetwork + DisGeNet datasets to perform GSEA/ORA
## kegg.do                logical, logical    Use the KEGG dataset to perform GSEA/ORA
## wp.do                  logical, logical    Use the WikiPathways dataset to perform GSEA/ORA
## reactome.do            logical, logical    Use the Reactome dataset to perform GSEA/ORA
## cm.do                  logical, logical    Use the CellMarker dataset to perform GSEA/ORA
## mesh.do                logical, logical    Use the MeSHDb dataset collection to perform GSEA/ORA (warning, this may consume an astronomical amount of RAM. gsea.do is forced to false)
## species                character           Name of the species analyzed (namely, 'human' or 'mouse')
## enrp.max               0<numeric<1         BH FDR-adjusted p-value cut-off to consider enriched terms as significant
## enr.min.genes          numeric+            Minimum number of significant genes to perform GSEA/ORA
## or.top.max             numeric+            Maximum number of significant genes to consider as a signature for ORA
## dotplot.maxterms       numeric+            Maximum number of enriched terms to plot in a dotplot (for readability)
## only.others            logical             When number of classes > 2, only perform "N_vs_Others" tests


DE.test <- function(exp.mat = NULL, annot.df = NULL, samples.colname = "Sample", batch.colnames = NULL, factorname = NULL, invert.levels = FALSE, adjp.max = 5E-02, lfc.min = 1, lfcShrink = TRUE, enrp.max = 1E-02, enr.min.genes = 10, or.top.max = 100, only.others = FALSE, outdir = getwd(), samples.dist.method = 'spearman', samples.hclust.method = 'ward.D', genes.dist.method = 'pearson', genes.hclust.method = 'ward.D', msigdb.do = c(TRUE, TRUE), do.do = c(TRUE, TRUE), go.do = c(TRUE, TRUE), kegg.do = c(TRUE, TRUE), wp.do = c(TRUE, TRUE), reactome.do = c(TRUE, TRUE), cm.do = c(TRUE, TRUE), mesh.do = c(FALSE, FALSE), non.redundent = TRUE, species = 'Homo sapiens', dotplot.maxterms = 50, my.seed = 1337, BPPARAM = BiocParallel::SerialParam()) {
  
  if (tolower(species) == 'homo sapiens') {
    Org <- 'org.Hs'
  } else if (tolower(species) == 'mus musculus') {
    Org <- 'org.Mm'
  } else stop("Only 'Homo sapiens' and 'Mus musculus' species are supported !")
  
  # message(cur.fac)
  
  ## CHECKS
  ## if exp.mat and annot.df have different size
  if (nrow(annot.df) != ncol(exp.mat)) stop("'exp.mat' and 'annot.df' do not have the same amount of samples!")
  ## if samples.colname does not exist
  if (!samples.colname %in% colnames(annot.df)) stop(paste0("Sample column '", samples.colname, "' not found !"))
  ## if any of batch.colnames does not exist
  if (!all(batch.colnames %in% colnames(annot.df))) stop("Some of batch column name(s) not found !")
  ## if cur.fac does not exist
  factn.check <- factorname %in% colnames(annot.df)
  if (!any(factn.check)) stop('None of the factor name(s) were found in the annotation table !')
  if (!all(factn.check)) {
    out.factn <- factorname[!factn.check]
    factorname <- factorname[factn.check]
    message(paste0('Some factor name(s) [', paste(out.factn, collapse = ', '), '] were not found in the annotation table, thus discarded'))
  }
  
  ## if samples are not identical
  exp.mat <- exp.mat[,order(colnames(exp.mat))]
  annot.df <- annot.df[order(annot.df[[samples.colname]]),]
  if (!all(colnames(exp.mat) == annot.df[[samples.colname]])) stop(paste0("Content of the sample column '", samples.colname, "' is not identical to 'exp.mat' colnames !"))
  
  library(DESeq2)
  library(SummarizedExperiment)
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
  
  ## Loading CellMarker db if needed
  if (any(cm.do)) {
    `%>%` <- dplyr::`%>%`
    cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>% tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% dplyr::select(cellMarker, geneID) %>% dplyr::mutate(geneID = strsplit(geneID, ', '))
  }
  
  ## Looping on factor names
  for (cur.fac in factorname) {
    
    ## Limiting annot.df to required columns
    cur.annot.df <- annot.df[, colnames(annot.df) %in% c(samples.colname, batch.colnames, cur.fac)]
    cur.exp.mat <- exp.mat
    cur.annot.df <- annot.df
    ## Adjusting the datasets if needed (handling NAs in annotation)
    na.check <- is.na(cur.annot.df[[cur.fac]])
    if (any(na.check)) {
      cur.exp.mat <- cur.exp.mat[,!na.check]
      cur.annot.df <- cur.annot.df[!na.check,]
    }
    
    ## Forcing a relevel (if some levels were lost)
    levtab <- as.data.frame(table(cur.annot.df[[cur.fac]]), stringsAsFactors = FALSE)
    levtab <- levtab[order(levtab$Var1),]
    levtab <- levtab[levtab$Freq > 0,]
    cur.annot.df[[cur.fac]] <- as.factor(as.character(cur.annot.df[[cur.fac]]))
    myref <- levels(cur.annot.df[[cur.fac]])[1]
    cur.annot.df[[cur.fac]] <- relevel(cur.annot.df[[cur.fac]], ref = myref)
    
    ## Filtering special characters
    levels(cur.annot.df[[cur.fac]]) <- gsub(pattern = "\\W", replacement = '.', x = levels(cur.annot.df[[cur.fac]]))
    
    ## Preparing the table list of all combinations
    mylevels <- levels(cur.annot.df[[cur.fac]])
    combn.res <- if(invert.levels) combn(mylevels, 2) else combn(rev(mylevels), 2)
    all.combz <- sapply(1:ncol(combn.res), function(x) {list(combn.res[1,x], combn.res[2,x])}, simplify = FALSE)
    names(all.combz) <- vapply(1:ncol(combn.res), function(x) { paste(combn.res[, x, drop = TRUE], collapse = '_vs_') }, 'a')
    if(length(mylevels) > 2) {
      XvsO.combz <- sapply(as.character(mylevels), function(x) { list(x, mylevels[!mylevels == x])}, simplify = FALSE)
      names(XvsO.combz) <- vapply(mylevels, function(x) { paste0(x, '_vs_Other') }, 'a')
      all.combz <- if(only.others) XvsO.combz else c(all.combz, XvsO.combz)
      rm(XvsO.combz)
    }
    
    ## Filtering comparisons with a class comprised by an unique sample
    for (mycomb in names(all.combz)) {
      class.len <- vapply(all.combz[[mycomb]], function(x) { length(which(cur.annot.df[[cur.fac]] %in% x))}, 1L)
      if(any(class.len == 1)) all.combz[[mycomb]] <- NULL
    }
    
    ## Creating design (factor then Batch, w/o intercept)
    my.textform <- paste(c('~0', unique(c(cur.fac, batch.colnames))), collapse = '+')
    my.design <- as.formula(my.textform)
    
    ## Creating the main output dir
    fac.dir <- paste0(outdir, '/Differential_analysis')
    factor.dir <- paste(c(fac.dir, paste0('adjp.', adjp.max, '_lfc.', lfc.min), my.textform), collapse = '/')
    dir.create(path = factor.dir, recursive = TRUE)
    
    ## Creating the DESeq2 object
    DE2obj <- DESeq2::DESeqDataSetFromMatrix(countData = cur.exp.mat, colData = cur.annot.df, design = my.design)
    rm(cur.exp.mat, cur.annot.df)
    
    ## Saving the DESeq object
    saveRDS(object = DE2obj, file = paste0(factor.dir, '/', cur.fac, '_rawcounts.RDS'), compress = 'bzip2')
    
    ## Normalizing by vst (for heatmap only)
    DE2obj.norm <- DESeq2::vst(object = DE2obj, blind = TRUE)
    norm.mat <- SummarizedExperiment::assay(DE2obj.norm)
    rm(DE2obj.norm)
    ### Removing batch effect if needed
    if (!is.null(batch.colnames)) {
      batch2.name <- if (length(batch.colnames) > 1) batch.colnames[2] else NULL
      norm.mat <- limma::removeBatchEffect(x = norm.mat, batch = SummarizedExperiment::colData(DE2obj)[[batch.colnames[1]]], batch2 = batch2.name)
    }

    ### PCAs
    for (p in c(cur.fac, batch.colnames)) {
      png(filename = paste0(factor.dir, '/PCA_vst_', p, '.png'), width = 1024, height = 1000)
      library(ggfortify)
      print(ggplot2::autoplot(prcomp(t(norm.mat)), data = as.data.frame(SummarizedExperiment::colData(DE2obj)), colour = p, size = 3))
      dev.off()
    }
    
    ## Performing the DE test
    htg.de.wald <- DESeq2::DESeq(DE2obj)
    
    ## Saving the Wald test DESeq object
    saveRDS(object = htg.de.wald, file = paste0(factor.dir, '/', cur.fac, '_wald.RDS'), compress = 'bzip2')
    
    ## LOOPING through pair combinations
    for (mycomb in names(all.combz)) {
      
      ## Creating output dir
      mycoef <- paste0(cur.fac, '~', mycomb)
      message(mycoef)
      
      de.dir <- paste(c(factor.dir, mycomb), collapse = '/')
      dir.create(path = de.dir, recursive = TRUE)
      
      ## Getting results table for current contrast
      mycontrast <- sapply(all.combz[[mycomb]], function(x) { paste0(cur.fac, x)}, simplify = FALSE)
      DEres <- DESeq2::results(htg.de.wald, contrast = mycontrast, listValues = c(1, -1/length(all.combz[[mycomb]][[2]])), independentFiltering = TRUE, alpha = adjp.max, pAdjustMethod = "BH", parallel = TRUE, BPPARAM = BPPARAM)
      
      ## Saving the test results object
      saveRDS(object = DEres, file = paste0(de.dir, '/', mycomb, '_results.RDS'), compress = 'bzip2')
      
      ## Shrinking l2fc
      if (lfcShrink) {
        suppressMessages(DEres <- DESeq2::lfcShrink(htg.de.wald, contrast = mycontrast, type = 'ashr', res = DEres))
        ## Saving the test reults object
        saveRDS(object = DEres, file = paste0(de.dir, '/', mycomb, '_results_lfcShrink.RDS'), compress = 'bzip2')
      }
      # rm(htg.de.wald)
      
      ## Histogram of P-values
      png(paste0(de.dir, '/', mycomb, '_phist.png'), width = 2048, height = 768)
      par <- par(mfrow = c(1, 2))
      hist(DEres$pvalue, col = "lightblue", main = paste0("Histogram of raw P-values (DESeq2)\n", mycoef), breaks = 100, xlim = c(0,1), xlab = "P-value")
      hist(DEres$padj, col = "lightblue", main = paste0("Histogram of BH-adjusted P-values (DESeq2)\n", mycoef), breaks = 100, xlim = c(0,1), xlab = "P-value")
      abline(v = adjp.max, col = 2, lty = 2)
      dev.off()
      
      ## MAplot
      png(paste0(de.dir, '/', mycomb, '_MA.png'), width = 1024, height = 768)
      DESeq2::plotMA(DEres, alpha = adjp.max, main = paste0("M-A Plot\n", mycoef), cex = 1)
      dev.off()
      
      ## Volcano plot
      deg.idx <- DEres$padj <= adjp.max & abs(DEres$log2FoldChange) >= lfc.min
      png(paste0(de.dir, '/', mycomb, '_volcano.png'), width = 1024, height = 768)
      plot(x = DEres$log2FoldChange, y = -log10(DEres$padj), xlab = "log2(Fold-Change)", ylab = "-log10(adjusted P-value)", col = ifelse(deg.idx, "red", "black"), main = paste0("Volcano plot\n", mycoef), pch = 20)
      grid()
      abline(h = -log10(adjp.max), lty = 2, col = 4)
      abline(v = lfc.min * c(-1, 1), lty = 2, col = 4)
      dev.off()
      ## Output table
      DEres.df <- cbind(Symbol = rownames(DEres), as.data.frame(DEres))
      DEres.df$cuts.in <- 0
      DEres.df$cuts.in[deg.idx] <- 1
      DEres.df <- DEres.df[order(DEres.df$padj, abs(DEres.df$log2FoldChange), decreasing = c(F, T)),]
      write.table(DEres.df, file = paste0(de.dir, '/', mycomb, '_results.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
      
      sig.genes <- as.character(DEres.df$Symbol[DEres.df$cuts.in == 1])
      
      ## Setting a color palette for the heatmaps
      myPalette <- c("royalblue3", "ivory", "orangered3")
      myRamp <- circlize::colorRamp2(c(-2, 0, 2), myPalette)
      
      if (length(sig.genes) > enr.min.genes) {
        
        cur.annot <- as.data.frame(SummarizedExperiment::colData(DE2obj)[,c(cur.fac, batch.colnames), drop = FALSE])
        cur.samples.idx <- cur.annot[[cur.fac]] %in% unlist(all.combz[[mycomb]])
        
        ## Heatmap
        ## data to plot
        plotDat <- norm.mat[rownames(norm.mat) %in% sig.genes, cur.samples.idx]
        # plotDat <- norm.mat[rownames(norm.mat) %in% sig.genes, SummarizedExperiment::colData(DE2obj)[[cur.fac]] %in% unlist(all.combz[[mycomb]])]
        z.mat <- (plotDat - rowMeans(plotDat)) / matrixStats::rowSds(plotDat)
        # Creating sample annotation
        # ha1 = ComplexHeatmap::HeatmapAnnotation(df = SummarizedExperiment::colData(DE2obj)[,c(cur.fac, batch.colnames), drop = FALSE][SummarizedExperiment::colData(DE2obj)[[cur.fac]] %in% unlist(all.combz[[mycomb]]),])
        ha1 = ComplexHeatmap::HeatmapAnnotation(df = cur.annot[cur.samples.idx,c(cur.fac, batch.colnames), drop = FALSE])
        
        ## Looping through requested clustering methods
        for (sdm in samples.dist.method) {
          for (shm in samples.hclust.method) {
            for (gdm in genes.dist.method) {
              for (ghm in genes.hclust.method) {
                ## Clustering samples
                hc.s <- hclust(amap::Dist(x = t(plotDat), method = sdm), method = shm)
                ## Clustering genes
                hc.g <- hclust(amap::Dist(x = plotDat, method = gdm), method = ghm)
                # Enhanced heatmap
                png(paste0(de.dir, '/', mycomb, '_sig.', nrow(z.mat), 'x', ncol(z.mat), '_', paste(c(gdm, ghm, sdm, shm), collapse = "_"), '.heatmap.png'), width = min(ncol(z.mat) * 15, 2000) + 200, height = min(length(sig.genes) * 10, 5000) + 300)
                suppressMessages(print(ComplexHeatmap::Heatmap(z.mat, name = "Normalized counts",
                                              # use my custom color palette
                                              col = myRamp,
                                              # do not show gene names
                                              show_row_name = TRUE,
                                              # do not clusterize samples
                                              cluster_columns = hc.s,
                                              cluster_rows = hc.g,
                                              # add a nice grey border to cells
                                              rect_gp = grid::gpar(col = "darkgrey", lwd=0.5),
                                              # add sample annotation
                                              top_annotation = ha1)))
                dev.off()
              }
            }
          }
        }
      }
      
      ## Shorter heatmap if more sig genes than requested "topN"
      if(length(sig.genes) > or.top.max) {
        
        ## Heatmap
        ## data to plot
        sig.genes <- as.character(DEres.df$Symbol[DEres.df$cuts.in == 1][1:or.top.max])
        plotDat <- norm.mat[rownames(norm.mat) %in% sig.genes, cur.samples.idx]
        # plotDat <- norm.mat[rownames(norm.mat) %in% sig.genes, SummarizedExperiment::colData(DE2obj)[[cur.fac]] %in% unlist(all.combz[[mycomb]])]
        z.mat <- (plotDat - rowMeans(plotDat)) / matrixStats::rowSds(plotDat)
        # Creating sample annotation
        # ha1 = ComplexHeatmap::HeatmapAnnotation(df = SummarizedExperiment::colData(DE2obj)[,c(cur.fac, batch.colnames)][annot.df[[cur.fac]] %in% unlist(all.combz[[mycomb]]),])
        ha1 = ComplexHeatmap::HeatmapAnnotation(df = cur.annot[cur.samples.idx,c(cur.fac, batch.colnames), drop = FALSE])
        ## Clustering samples
        for (sdm in samples.dist.method) {
          for (shm in samples.hclust.method) {
            for (gdm in genes.dist.method) {
              for (ghm in genes.hclust.method) {
                hc.s <- hclust(amap::Dist(x = t(plotDat), method = sdm), method = shm)
                ## Clustering genes
                hc.g <- hclust(amap::Dist(x = plotDat, method = gdm), method = ghm)
                # Enhanced heatmap
                png(paste0(de.dir, '/', mycomb, '_sig.TOP', nrow(z.mat), 'x', ncol(z.mat), '_', paste(c(gdm, ghm, sdm, shm), collapse = "_"), '.heatmap.png'), width = min(ncol(z.mat) * 15, 2000) + 200, height = min(length(sig.genes) * 10, 5000) + 300)
                suppressMessages(print(ComplexHeatmap::Heatmap(z.mat, name = "Normalized counts",
                                              # use my custom color palette
                                              col = myRamp,
                                              # do not show gene names
                                              show_row_name = TRUE,
                                              # do not clusterize samples
                                              cluster_columns = hc.s,
                                              cluster_rows = hc.g,
                                              # add a nice grey border to cells
                                              rect_gp = grid::gpar(col = "darkgrey", lwd=0.5),
                                              # add sample annotation
                                              top_annotation = ha1)))
                dev.off()
              }
            }
          }
        }
      }
      
      ## Functional enrichment
      if (any(msigdb.do, kegg.do, do.do, go.do, cm.do, wp.do, reactome.do, mesh.do) & length(which(deg.idx)) >= enr.min.genes & !is.null(species)) {
        
        enr.inputs <- table2enr(deseq2.res.data = DEres.df, species = species, geneid.colname = 'Symbol', geneid.type = 'SYMBOL', value.colname = 'log2FoldChange', topN.max = or.top.max, topN.order.colname = 'padj', topN.order.decreasing = FALSE, topN.cutoff = enrp.max, topN.keep.operator = '<')
        
        ## MSIGDB
        if (any(msigdb.do)) {
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
              if (!is(my.gsea.res, class2 = 'try-error')) gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = mycomb)
            }
            ### ORA
            if(msigdb.do[2]) {
              my.ora.res <- try(ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = 'clusterProfiler::enricher', t2g = my.t2g, t2g.name = msc, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
              ## Generate plots / outputs
              if (!is(my.ora.res, class2 = 'try-error')) ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = mycomb, geneList = enr.inputs$gsea.genevec)
            }
          }
        }
        
        ## GO (gene ontology)
        if (any(go.do)) {
          ### GSEA
          if(go.do[1]) {
            func.name <- 'clusterProfiler::gseGO'
            for (x in c('BP', 'CC', 'MF')) {
              my.gsea.res <- try(gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enrp.max, minGSSize = enr.min.genes, OrgDb = get(paste0(msigdbr2org(species), '.db')), ont = x))
              if (!is(my.gsea.res, class2 = 'try-error')) {
                my.gsea.res@setType <- paste(c(my.gsea.res@setType, x), collapse = '_')
                gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = mycomb)
                ## Simplify
                if(nrow(my.gsea.res) > 1) {
                  my.gsea.res@setType <- x
                  my.gsea.res <- enrichplot::pairwise_termsim(my.gsea.res)
                  my.gsea.res.simp <- clusterProfiler::simplify(my.gsea.res, cutoff = 0.7, by = "p.adjust", select_fun = min)
                  if(nrow(my.gsea.res.simp) < nrow(my.gsea.res)) {
                    my.gsea.res.simp@setType <- paste(c(func.name, x, 'simplified'), collapse = '_')
                    gsea.output(gseaResult = my.gsea.res.simp, out.dir = de.dir, comp.name = mycomb)
                  }
                }
              }
            }
          }
          ### ORA
          if(go.do[2]) {
            func.name <- 'clusterProfiler::enrichGO'
            for (x in c('BP', 'CC', 'MF')) {
              my.ora.res <- ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes, OrgDb = get(paste0(msigdbr2org(species), '.db')), ont = x)
              if (!is(my.ora.res, class2 = 'try-error')) {
                my.ora.res@ontology <- paste(c(my.ora.res@ontology, x), collapse = '_')
                ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = mycomb, geneList = enr.inputs$gsea.genevec)
                ## Simplify
                if(nrow(my.ora.res) > 1) {
                  my.ora.res@ontology <- x
                  my.ora.res <- enrichplot::pairwise_termsim(my.ora.res)
                  my.ora.res.simp <- clusterProfiler::simplify(my.ora.res, cutoff = 0.7, by = "p.adjust", select_fun = min)
                  if(nrow(my.ora.res.simp) < nrow(my.ora.res)) {
                    my.ora.res.simp@ontology <- paste(c(func.name, x, 'simplified'), collapse = '_')
                    ora.output(enrichResult = my.ora.res.simp, out.dir = de.dir, comp.name = mycomb, geneList = enr.inputs$gsea.genevec)
                  }
                }
              }
            }
          }
        }
        
        ## DO (disease ontology)
        if (any(do.do)) {
          ### GSEA
          if (do.do[1]) {
            for (x in c('DOSE::gseDO', 'DOSE::gseNCG', 'DOSE::gseDGN')) {
              my.gsea.res <- try(gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
              if (!is(my.gsea.res, class2 = 'try-error')) gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = mycomb)
            }
            ### ORA
            if(do.do[2]) {
              for (x in c('DOSE::enrichDO', 'DOSE::enrichNCG', 'DOSE::enrichDGN')) {
                my.ora.res <- try(ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
                if (!is(my.ora.res, class2 = 'try-error')) ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = mycomb, geneList = enr.inputs$gsea.genevec)
              }
            }
          }
        }
        
        ## KEGG/MKEGG
        ### NOTE1 : It's the same way to call the 'gsea.run' / 'ora.run' as it is for 'DO', 'NCG' or 'DGN', but here it's compatible with many more species than homo sapiens.
        ### NOTE2 : for this case, additional KEGG pathway plots will be generated.
        ### NOTE3 : for this case, an internet connexion is required to query the KEGG website.
        if (any(kegg.do)) {
          ### GSEA
          if (kegg.do[1]) {
            for (x in c('clusterProfiler::gseKEGG', 'clusterProfiler::gseMKEGG')) {
              my.gsea.res <- try(gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
              if (!is(my.gsea.res, class2 = 'try-error')) gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = mycomb)
            }
            ### ORA
            if(kegg.do[2]) {
              for (x in c('clusterProfiler::enrichKEGG', 'clusterProfiler::enrichMKEGG')) {
                my.ora.res <- try(ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = x, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
                if (!is(my.ora.res, class2 = 'try-error')) ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = mycomb, geneList = enr.inputs$gsea.genevec)
              }
            }
          }
        }
        
        ## WP (wikipathways)
        if(any(wp.do)) { 
          if(wp.do[1]) {
            ### GSEA
            func.name <- 'clusterProfiler::gseWP'
            my.gsea.res <- try(gsea.run(geneList = enr.inputs$gsea.genevec, organism = species, func.name = func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
            if (!is(my.gsea.res, class2 = 'try-error')) gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = mycomb)
          }
          ### ORA
          if(wp.do[2]) {
            func.name <- 'clusterProfiler::enrichWP'
            my.ora.res <- ora.run(gene = enr.inputs$ora.genevec, organism = species, func.name = func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes)
            if (!is(my.ora.res, class2 = 'try-error')) ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = mycomb, geneList = enr.inputs$gsea.genevec)
          }
        }
        
        ## REACTOME
        if (any(reactome.do)) {
          reactome.org <- tolower(convert_species_name(OrgDb = get(paste0(msigdbr2org(species = species), '.db'))))
          if(reactome.do[1]) {
            ### GSEA
            func.name <- 'ReactomePA::gsePathway'
            my.gsea.res <- gsea.run(geneList = enr.inputs$gsea.genevec, organism = reactome.org, func.name = func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enrp.max, minGSSize = enr.min.genes)
            my.gsea.res@setType <- paste0(func.name, '_Reactome')
            gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = mycomb)
          }
          if(reactome.do[2]) {
            ### ORA
            func.name <- 'ReactomePA::enrichPathway'
            my.ora.res <- ora.run(gene = enr.inputs$ora.genevec, organism = reactome.org, func.name = func.name, t2g = NULL, t2g.name = NULL, gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes)
            my.ora.res@ontology <- paste0(func.name, '_Reactome')
            ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = mycomb, geneList = enr.inputs$gsea.genevec)
          }
        }
        
        ## CELLMARKER
        ### Assess cell types from an online table
        ### NOTE : It's the same way to call the 'gsea.run' / 'ora.run' functions as for MSIGDB, but with a single bank (so, no loop).
        if (any(cm.do)) {
          ### GSEA
          if (cm.do[1]) {
            my.gsea.res <- try(gsea.run(geneList = enr.inputs$gsea.genevec, species = species, func.name = 'clusterProfiler::GSEA', t2g = cell_markers, t2g.name = 'CellMarker', gene2Symbol = enr.inputs$gene2Symbol, seed = my.seed, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
            if (!is(my.gsea.res, class2 = 'try-error')) gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = mycomb)
          }
          #### ORA
          if (cm.do[2]) {
            my.ora.res <- try(ora.run(gene = enr.inputs$ora.genevec, species = species, func.name = 'clusterProfiler::enricher', t2g = cell_markers, t2g.name = 'CellMarkers', gene2Symbol = enr.inputs$gene2Symbol, pvalueCutoff = enrp.max, minGSSize = enr.min.genes))
            if (!is(my.ora.res, class2 = 'try-error')) ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = mycomb, geneList = enr.inputs$gsea.genevec)
          }
        }
        
        ## MESH (WARNING : MEMORY OGRE AND SLOW !! Big DBs, 3 sources, 16 categories ! 64 GB of RAM required for most bases !
        ### Requires additional parameters :
        ### . 'MeSHDb' : character ; name of a MeSH [NO : AUTO FROM SPECIES NAME]
        ### . 'database' : character ; MeSH source type (can be 'gendoo' = text-mining, 'gene2pubmed' = manual curation by NCBI team, 'RBBH' = sequence homology with BLASTP search @ E-value < 1E-50)
        ### . 'category' : character ; name of a MeSH category sub-db (namely 'A', 'B', 'C', 'D', 'G').
        ### NOTE : see https://yulab-smu.top/biomedical-knowledge-mining-book/meshes-semantic-similarity.html
        if (any(mesh.do)) {
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
                    ora.output(enrichResult = my.ora.res, out.dir = de.dir, comp.name = mycomb, geneList = enr.inputs$gsea.genevec)
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
                    gsea.output(gseaResult = my.gsea.res, out.dir = de.dir, comp.name = mycomb)
                  }
                } else message(paste0("Unsupported MeSH database '", y, "'. Expecting one of : '", paste(mesh.dbs, collapse = "', '"), "'."))
              }
            }
          }
        }
      }
    }
  }
}

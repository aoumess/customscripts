## This function performs removal of ambient RNA expression for a single cell counts matrix using SoupX. Usage of the raw (ie, without discarding the empty barcodes) is recommended. Without, lower efficiency is expected. TESTED WITH SoupX v1.6.2
### scmat_filt      [num matrix]        Count matrix corresponding to barcodes kept as "true" cells (by DropletUtils::emptyDrops). [REQUIRED]
### scmat_raw       [num matrix]        Count matrix corresponding to all barcodes. [OPTIONAL, RECOMMENDED]
### Other parameters : see SoupX::autoEstCont
RunSoupX <- function(scmat_filt = NULL, scmat_raw = NULL, maxMarkers = 100, forceAccept = FALSE, contaminationRange = c(.01,.8)) {
  ## If no raw matrix
  if (is.null(scmat_raw)) {
    spChanRaw <- SoupX::SoupChannel(tod = scmat_filt, toc = scmat_filt, calcSoupProfile = FALSE)
    sc_rowsum <- Matrix::rowSums(scmat_filt)
    spProf <- data.frame(row.names = rownames(scmat_filt), est = sc_rowsum/sum(scmat_filt), counts = sc_rowsum)
    spChan <- SoupX::setSoupProfile(spChanRaw, spProf)
  } else {
    spChan <- SoupX::SoupChannel(tod = scmat_raw, toc = scmat_filt, calcSoupProfile = TRUE)
  }
  ## Quick clustering needed
  spClust <- scran::quickCluster(scmat_filt, method = "igraph")
  ## Adding clusters to the SoupChannel object
  spChan <- SoupX::setClusters(sc = spChan, clusters = spClust)
  ## Estimating soup
  sX <- SoupX::autoEstCont(sc = spChan, doPlot = FALSE, maxMarkers = maxMarkers, forceAccept = forceAccept, contaminationRange = contaminationRange)
  ## Removing soup (adjusting counts)
  scmat_soupx <- SoupX::adjustCounts(sX, method = 'subtraction', roundToInt = TRUE, tol = .001, pCut = .01)
  return(scmat_soupx)
}

## Example 1
# scmat_filt <- readRDS('/shared/projects/form_2022_32/SingleCellRNASeq/Platdum_post_droplet_filtered.RDS')
# scmat_raw <- readRDS('/shared/projects/form_2022_32/SingleCellRNASeq/Platdum_post_raw.RDS')
# soupx_res <- RunSoupX(scmat_filt = scmat_filt, scmat_raw = scmat_raw)
# saveRDS(object = soupx_res, file = '/shared/projects/form_2022_32/SingleCellRNASeq/Platdum_SoupXadjusted.RDS', compress = "bzip2")
# soupx_res_basic <- RunSoupX(scmat_filt = scmat_filt, scmat_raw = NULL)
# saveRDS(object = soupx_res_basic, file = '/shared/projects/form_2022_32/SingleCellRNASeq/Platdum_SoupXadjusted_basic.RDS', compress = "bzip2")

## Example 2
### Loading data
# load('/home/job/WORKSPACE/COLLABORATION/With_RemyJELIN/2303_GE_QC_NON-NORMALIZED.rda')
# scmat_filt <- sobj@assays$RNA@counts
# rm(sobj)
### We are expecting a very high contamination (emptyDrops estimates 48K cells from 10X chromium data were 20K cells were loaded), so we increase the amount of markers, as well as the contamination range, and force the output without getting into error. 
# scmat_soupx <- RunSoupX(scmat_filt = scmat_filt, maxMarkers = 500, forceAccept = TRUE, contaminationRange = c(.01,.95))
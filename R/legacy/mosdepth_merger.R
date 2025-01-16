workdir <- "/mnt/data_cigogne/job/B18_JBMI/Batch_02/09_QC/"

mdbed.list <- list.files(path = workdir, pattern = "\\.bed$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE, include.dirs = FALSE)

require(foreach)


depth.all <- foreach(mdbed = mdbed.list, .combine = "cbind") %do% {
  mdb.content <- read.table(mdbed, header = FALSE, sep = "\t", as.is = TRUE)
  if (mdbed == mdbed.list[1]) return(mdb.content) else return(mdb.content[,ncol(mdb.content)])
}

samplenames <- vapply(basename(mdbed.list), function(x) { paste(unlist(strsplit(x, split = "_"))[c(1,2)], collapse = "_")}, 'a')
colnames(depth.all) <- c("Chr", "Start", "End", "Target", samplenames)
write.table(depth.all, file = paste0(workdir, "/mosdepth_merged.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

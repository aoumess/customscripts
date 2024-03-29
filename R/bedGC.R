bed2gc <- function(bed = NULL, header = FALSE, genome.fasta.file = NULL, species = "human", genome.build = "hg19", nt.add = 0, nthread = 1) {

  ## Checks
  if(is.null(bed)) stop("Please provide a valid path to a BED file !")
  if(is.null(genome.fasta.file)) stop("Please provide a valid path to a genome FASTA file !")
  if(!is.numeric(nt.add)) stop("'nt.add' must be a numeric vector (actually an integer vector) !")
  if(any(nt.add < 0)) stop("'nt.add' must be a POSITIVE numeric vector !")
  if(any(is.infinite(nt.add))) stop("'nt.add' must be a FINITE positive numeric vector !")
  if(!(tolower(species) %in% c("human", "mouse"))) stop("'species' must be 'human' or 'mouse' !")
  if(!file.exists(bed)) stop(paste0("Could not open BED file ", bed, " !"))
  if(!file.exists(genome.fasta.file)) stop(paste0("Could not open genome FASTA file ", genome.fasta.file, " !"))
  require(parallel)
  logicores <- parallel::detectCores(logical = TRUE)
  if(nthread > logicores) print(paste0("WARNING ! ", nthread, "threads requested but only ", logicores, " logical threads available !"))

  ## Adjusting nt.add
  nt.add <- unique(nt.add)
  
  ## Adjusting threads
  if (length(nt.add) < nthread) nthread <- length(nt.add)

  ## Reading BED
  print("Loading BED file ...")
  require(maelstrom)
  my.bed <- maelstrom::read.table.fast(file = bed, header = header, sep = "\t")
  my.bed$start <- my.bed$start - 1  ## Tweak for bedtools nuc
  bed.ncol <- ncol(my.bed)
  if (bed.ncol == 3) colnames(my.bed) <- c("chr", "start", "end") else colnames(my.bed)[1:4] <- c("chr", "start", "end", "name")

  ## Preparing the BED data
  print("Formatting BED file ...")
  require(bedr)
  my.bed.sorted <- bedr(engine = "bedtools", input = list(i = my.bed), method = "sort", params = "", check.zero.based = FALSE, check.sort = FALSE, check.merge = FALSE)
  
  ## Preparing MT
  print(paste0("Loading cluster (", nthread, ") ..."))
  require(foreach)
  require(doParallel)
  cl <- makeCluster(spec = nthread, type = "FORK")
  registerDoParallel(cl)

  print("Computing GC tracks ...")
  myXgc <- foreach(nta = nt.add, .combine = "cbind") %dopar% {
    print(nta)
    if (nta > 0)  my.bed.grown <- grow.region(x = my.bed.sorted, n.add = nta, species = species, build = genome.build, check.zero.based = FALSE, check.sort = FALSE, check.merge = FALSE, verbose = FALSE) else my.bed.grown <- my.bed.sorted
    my.gc <- bedr(engine = "bedtools", input = list(bed = my.bed.grown), method = "nuc", params = paste0("-fi ", genome.fasta.file), check.merge = FALSE, check.sort = FALSE, verbose = FALSE)
    mygc <- as.numeric(my.gc[-1,(ncol(my.bed.grown)+2)])
    return(mygc)
  }

  print("Stopping cluster ...")
  stopCluster(cl)

  mypos <-cbind(my.bed, myXgc)
  colnames(mypos)[(bed.ncol+1):ncol(mypos)] <- paste0(nt.add, "bp")
  

  print("Writing GC data ...")
  maelstrom::write.table.fast(x = mypos, file = paste0(bed, ".gc"), header = TRUE, sep = "\t")
  print("Done.")
}

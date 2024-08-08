work_dir <- '/home/job/WORKSPACE/PROJECTS/B22011_VAMA_01/customplot/'
upin <- 'broad_consensus_peaks.mRp.clN.boolean.intersect.txt'
upin <- 'narrow_consensus_peaks.mRp.clN.boolean.intersect.txt'
intersect_min <- 0
pdf.factor <- 2

## Load data
infile <- paste0(work_dir, '/', upin)
updf <- read.table(file = infile, header = FALSE, sep = '\t', as.is = TRUE)
updat <- setNames(object = updf[,2], nm = updf[,1])

## Create upset df
upx <- UpSetR::fromExpression(updat[updat >= intersect_min])

sets <- list(
  '1' = c('HEALTHY_CTRL', 'HEALTHY_IL6', 'HEALTHY_calpro', 'HEALTHY_IL6calpro')
  , '2' = c('HEALTHY_CTRL', 'HEALTHY_IL6')
  , '3' = c('MYELO_CTRL', 'MYELO_IL6', 'MYELO_calpro', 'MYELO_IL6calpro')
  , '4' = c('MYELO_CTRL', 'MYELO_IL6')
)

outfile <- sub(pattern = '\\.txt$', replacement = '.pdf', x = infile)
pdf(file = outfile, width = 29.7/cm(pdf.factor), height = 21/cm(pdf.factor))
for (setx in sets) print(UpSetR::upset(upx, order.by = 'freq', keep.order = TRUE, nsets = 100, nintersects = NA, sets = setx, number.angles = 45))
dev.off()




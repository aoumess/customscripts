## Load one or multiple GeoMX DSP xlsx output file (at the background correction step).
## This file should contain at least the 'SegmentProperties' panel (which contain QC metrics at the ROI/AOI level), the 'TargetCountMatrix' panel (which contains background-corrected raw counts) and the 'Dataset summary' panel (which exhibits which data QC/modifications were performed in DSP).
## All 3 tables are extracted for each provided path
## . geomx_files      [vec(char)]       Path(s) to the DSP xlsx file(s)
## . split_qcflags    [logical]         If QC multiple flags are found, split them to multiple new columns
## . qcflags_colname  [char]            Name of the column to find QC flags
## . batch_colname    [char]            How to name the batch column (created when multiple geomx_files are provided)
## . add_aoilabel     [logical]         When repeated ROI labels are found for a same slide, create the AOILabel column
geomx_xlsx_parser <- function(geomx_files = NULL, split_qcflags = TRUE, qcflags_colname = 'QCFlags', batch_colname = 'Batch', add_aoilabel = TRUE) {
  ## Check params
  if (is.null(geomx_files)) stop('At least one input files is required !')
  if (!is.logical(split_qcflags)) stop('The "split_qcflags" parameter should be logical !')
  if (!is.character(batch_colname)) stop('The output batch column name should be a character !')
  ## Check files
  if (!all(vapply(geomx_files, file.exists, TRUE))) stop('At least one of the provided GeoMX input file(s) do(es) not exist !')
  
  panel.names <- c('SegmentProperties', 'TargetCountMatrix', 'Dataset summary')
  
  ## SegmentProperties
  ### Reading the panel
  sp_list <- lapply(geomx_files, function(xl) { as.data.frame(readxl::read_excel(path = xl, sheet = panel.names[1], na = c('', 'na', 'NA'), trim_ws = TRUE, progress = FALSE)) })
  ## Check repeted ROI levels
  if (add_aoilabel) {
    for (xl in seq_along(geomx_files)) {
      slidenames <- unique(sp_list[[xl]][['SlideName']])
      repcheck <- any(sapply(unique(sp_list[[xl]][['SlideName']]), function(x) { any(duplicated(sp_list[[xl]][sp_list[[xl]][['SlideName']] == x,'ROILabel'])) }))
      if(!repcheck) {
        message('No replicated ROI found in ROILabel, so no addition of AOILabel.')
      } else {
        sp_list[[xl]][['AOILabel']] <- ''
        for (sn in slidenames) {
          snidx <- which(sp_list[[xl]][['SlideName']] == sn)
          sp_list[[xl]][['AOILabel']][snidx] <- sprintf("%03i", seq.int(1,length(snidx)))
        }
      }
    }
  }
  ### Adding the batch number
  for (xl in seq_along(geomx_files)) sp_list[[xl]][[batch_colname]] <- xl
  ### Merging (WARNING, the files may have different column amount and names !)
  sp_df <- Reduce(f = function(x,y) merge(x, y, all = TRUE), sp_list)
  ## Re-sorting
  if ('AOILabel' %in% colnames(sp_df)) sort_col <- 'AOILabel' else 'ROILabel'
  sp_df <- sp_df[order(sp_df[['SlideName']], sp_df[[sort_col]]),]
  rm(sp_list)
  ### Splitting QC flags when requested
  if (split_qcflags) {
    flags.list <- lapply(sp_df[[qcflags_colname]], function(x) { unlist(strsplit(x = x, split = ',')) })
    cur.flags <- sort(unique(unlist(flags.list)))
    if (length(cur.flags) == 0) {
      message('No QC flag to merge !')
    } else {
      for (fl in cur.flags) {
        fl2 <- paste(c('QC', unlist(strsplit(x = fl, split = '\\s+', fixed = FALSE))), collapse = '_')
        sp_df[[fl2]] <- vapply(flags.list, function(x) fl %in% x, TRUE)
        rm(fl2)
      }
      rm(flags.list)
    }
  }
  
  ## TargetCountMatrix
  ### Reading the panel
  tc_list <- lapply(geomx_files, function(xl) { as.data.frame(readxl::read_excel(path = xl, sheet = panel.names[2], na = c('', 'na', 'NA'), trim_ws = TRUE, progress = FALSE)) })
  ### Merging 
  tc_df <- Reduce(function(x,y) { y[['TargetName']] <- NULL; cbind(x, y) }, tc_list)
  rm(tc_list)
  tc_mat <- as.matrix(x = tc_df[,-1])
  dimnames(tc_mat) <- list(tc_df$TargetName, colnames(tc_df)[-1])
  rm(tc_df)
  
  ## Dataset summary
  ### Read the panel
  ds_list <- lapply(geomx_files, function(xl) { as.data.frame(readxl::read_excel(path = xl, sheet = panel.names[3], na = c('', 'na', 'NA'), trim_ws = TRUE, progress = FALSE, col_names = c('Key', 'Value'))) })
  names(ds_list) <- geomx_files
  
  ## Return list
  out.list <- list(sp_df, tc_mat, ds_list)
  names(out.list) <- panel.names
  rm(sp_df, tc_mat, ds_list)
  return(out.list)
  
}

## Test
# geomx_files <- c('/home/job/WORKSPACE/PROJECTS/B23025_MIBO_01_RT09722_GeoMx/DATA/Run_01_Test/RT09722_Run_01-Test_BGsub_BC_geomean_20230414.xlsx', '/home/job/WORKSPACE/PROJECTS/B23025_MIBO_01_RT09722_GeoMx/DATA/Run_02/RT09722_Run_02_BGsub_BC_geomean_20230821.xlsx')
# split_qcflags <- TRUE
# qcflags_colname <- 'QCFlags'
# batch_colname <- 'Batch'
# add_aoilabel <- TRUE

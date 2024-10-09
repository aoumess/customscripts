## Convert SVG to PNG ====
.svg_convert <- function(svg_files = NULL, format = 'png', compress = TRUE, ...) {
  valid_formats <- c('eps', 'png', 'pdf', 'ps', 'svg')
  if (!format %in% valid_formats) stop(paste0('Unsupported format ! Supported formats are : "', paste(valid_formats, collapse = '", "'), '"'))
  for (s in svg_files) {
    do.call(what = eval(parse(text = paste0('rsvg::rsvg_', format))), args = list(svg = s, file = sub(pattern = '\\.svg$', replacement = paste0('.', tolower(format)), x = s, ignore.case = TRUE), ...))
    if (compress) R.utils::gzip(s, remove = TRUE, overwrite = TRUE)
  }
}

## SVG2PNG+GZ ====
## Device manager to convert an open 'svg' device to raster (PNG by default) and compress (gzip) the original SVG, then close the device
svg_off <- function(format = 'png', compress = TRUE, ...) {
  ## Checks
  valid_formats <- c('eps', 'png', 'pdf', 'ps', 'svg')
  if (!format %in% valid_formats) stop(paste0('Unsupported format ! Supported formats are : "', paste(valid_formats, collapse = '", "'), '"'))
  ## Get current device
  cur_dev <- .Device
  ## Check if device is actually a SVG
  if (!cur_dev == 'svg') stop(paste0("Latest open device is not 'svg' : ", as.character(cur_dev)))
  ## Get filename
  svg_file <- attr(x = cur_dev, which = 'filepath', exact = TRUE)
  ## Get device dimensions
  dev_w <- dev.size()[1]*96
  dev_h <- dev.size()[2]*96
  ## Close device
  dev.off()
  ## Convert to raster
  .svg_convert(svg_files = svg_file, format = format, compress = compress, width = dev_w, height = dev_h)
}

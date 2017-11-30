#' Download raw files from GEO for a given study
#'
#' @param x a GEO Series id (GSE)
#' @param y the path to save the downloaded files. By default this value is set to the working directory
#' @return  downloadRaw creates a directory in the given path with the GSE name and saves there the compressed RAW files.
#' @examples
#'
#' Download the raw files of GSE11761 study and save them in a directory named GSE11761 in the working directory
#'
#' downloadRaw("GSE11761",getwd())
#'
#' @export
#' @importFrom GEOquery getGEOSuppFiles

downloadRaw <- function(x,y=getwd()) {
  if (missing(x))
    stop("Need to specify a GEO Series id, i.e 'GSE10026'")
  if ( (missing(y)) | (!file.exists(y)) )
    stop("Need to specify a valid path, i.e getwd()")
  if (!grepl("GSE[0-9]+",x))
    stop("x must be a GEO Series id, i.e 'GSE10026'")
  raw <- NULL
  a <- TRUE
  while(a) {
    raw <- try(GEOquery::getGEOSuppFiles(x, baseDir =y))
    if(!inherits(raw,"try-error")) a <- FALSE
    else cat("NCBI was busy, trying again\n")
  }

}

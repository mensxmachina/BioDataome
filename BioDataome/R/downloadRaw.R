#' Download raw CEL files from GEO for a given study
#'
#' @param x a GEO Series id (GSE)
#' @param y the path to save the downloaded files
#' @return a directory with the GSE with the compressed RAW files
#' @examples
#' downloadRaw("GSE11761")
#' @export
#' @importFrom GEOquery getGEOSuppFiles

downloadRaw <- function(x,y) {
  raw <- NULL
  a <- TRUE
  while(a) {
    raw <- try(GEOquery::getGEOSuppFiles(x, baseDir =y))
    if(!inherits(raw,"try-error")) a <- FALSE
    else cat("NCBI was busy, trying again\n")
  }

}

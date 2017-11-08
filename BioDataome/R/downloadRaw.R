#' Download raw CEL files from GEO for a given study
#'
#' @param x a GEO Series (GSE)
#' @return a directory with the GSE with the compressed RAW files
#' @examples
#' downloadRaw("GSE11761")
#' @export
#' @importFrom GEOquery
downloadRaw <- function(x) {
  raw <- NULL
  a <- TRUE
  while(a) {
    raw <- try(GEOquery::getGEOSuppFiles(x))
    if(!inherits(raw,"try-error")) a <- FALSE
    else cat("NCBI was busy, trying again\n")
  }

}

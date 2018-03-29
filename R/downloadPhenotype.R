#' Download series matrices from GEO for a given study
#'
#' This function downlads all series matrices related to a given GEO Series (GSE)
#' and saves them in a list.
#' The same GSE study (i.e. GSE11761) may be related  to more than one platforms (i.e GPL570 and GPL1261).
#' The length of the output list is the number of the related platforms.
#'
#' @param x a GEO Series id (GSE)
#' @return a list of series matrices related to the given study
#' @examples \dontrun{
#' downloadPhenotype("GSE11761") }
#' @export
#' @importFrom GEOquery getGEO

downloadPhenotype <- function(x) {
   if (missing(x))
    stop("Need to specify a GEO Series id, i.e 'GSE10026'")
  if (!grepl("GSE[0-9]+",x))
    stop("x must be a GEO Series id, i.e 'GSE10026'")
  phenotype <- NULL
  a <- TRUE
  while(a) {
    phenotype <- try(GEOquery::getGEO(x, GSEMatrix=TRUE,getGPL=FALSE))
    if(!inherits(phenotype,"try-error")) a <- FALSE
    else cat("NCBI was busy, trying again\n")
  }
  return(phenotype)
}

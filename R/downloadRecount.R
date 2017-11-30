#' Download gene-level RangedSummarizedExperiment data from Recount
#'
#' This function downloads the RangedSummarizedExperiment object with the data summarized at the gene level
#' from Recount (https://jhubiostatistics.shinyapps.io/recount/)
#'
#' @param x a recount dataset ID
#' @param y the destination path for the downloaded RangedSummarizedExperiment object
#' @return RangedSummarizedExperiment object for the given study
#' @examples
#' downloadRecount("SRP032775",getwd())
#' @export
#' @importFrom recount download_study

downloadRecount <- function(x) {
  rse_gene <- NULL
  a <- TRUE
  while(a) {
    rse_gene <- try(recount::download_study(x))
    if(!inherits(rse_gene,"try-error")) a <- FALSE
    else cat("Recount was busy, trying again\n")
  }
  return(rse_gene)
}

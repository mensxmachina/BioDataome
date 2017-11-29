#' Column wise comparison of two datasets
#' This function finds how many samples are shared between two datasets. The datasets are in the form
#' variables (probes) x samples. The number of variables (probes) in both datasets should be the same
#'
#' @param d1 a numeric matrix of a dataset
#' @param d2 a numeric matrix of a dataset
#' @return the number of equal samples
#' @examples
#' Let us assume we want to compare two normalized gene expression datasets from the same platform
#' d1<-get(load(url("http://dataome.mensxmachina.org/data/Homo%20sapiens/GPL570/GSE86013.Rda")))
#' d2<-get(load(url("http://dataome.mensxmachina.org/data/Homo%20sapiens/GPL570/GSE86015.Rda")))
#' commons<-compareDsets(d1,d2)
#' @export
#' @importFrom Rfast mat.mat

compareDsets<-function(d1,d2){

  if (missing(d1))
    stop("Need to specify a numeric matrix'")

  if (missing(d2))
    stop("Need to specify a numeric matrix'")

  if (!is.numeric(d1))
    stop("input matrix must be numeric")

  if (!is.numeric(d2))
  stop("input matrix must be numeric")

  if (nrow(d1)!=nrow(d2))
    stop("input matrices must have the same number of rows")

  a<-Rfast::mat.mat(d1,d2)
  b<-which(a==1)
  return(length(b))
}

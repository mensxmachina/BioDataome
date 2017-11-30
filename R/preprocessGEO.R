#' Preprocess CEL files with SCAN
#'
#' This function calls the SCAN method as described in Piccolo SR, Sun Y, Campbell JD, Lenburg ME,
#' Bild AH and Johnson WE (2012). A single-sample microarray normalization method to facilitate personalized-medicine workflows.
#' Genomics, 100(6), pp. 337-344.
#'
#' @param x the path where the CEL files are stored
#' @param y the number of cores to run in parallel
#' @return a matrix of dimensions: probes x samples with the normalized expression values
#' @examples
#' Assuming that CEL files are located in working directory
#' preprocessGEO(getwd(),3)
#' @export
#' @importFrom SCAN.UPC SCAN
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom Biobase exprs

preprocessGEO<-function(x,y){
  #get current directory to return after normalization
  currentPath<-getwd()
  setwd(x)
  cl <- parallel::makeCluster(y)
  doParallel::registerDoParallel(cl)
  normalized =SCAN.UPC::SCAN( "*.CEL.gz")
  parallel::stopCluster(cl)
  normalized<-Biobase::exprs(normalized)
  setwd(currentPath)
  return(normalized)
}



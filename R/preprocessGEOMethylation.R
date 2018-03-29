#' Preprocess IDAT files from Illumina HumanMethylation450 BeadChip
#'
#' This function utilizes minfi Package to convert data into methylation measurements.
#' Example: Assuming there is a directory named GSE78279 in the working
#' directory where idat files are stored and "GSM2071074_8655685078_R03C02"
#' and "GSM2071074_8655685078_R03C02" are the file names for the idat files,
#' then x should be:
#'
#' @param x a character array with the paths to the idat files
#' @return a matrix of dimensions: probes x samples with the normalized methylation values
#' @examples \dontrun{
#' x<-c("GSM2071074_8655685078_R03C02","GSM2071075_8655685078_R04C02")
#' x<-file.path(getwd(),"GSE78279",x)
#' dataNorm<-preprocessGEOMethylation(x) }
#' @export
#' @importFrom minfi read.metharray mapToGenome getBeta

preprocessGEOMethylation<-function(x){
  #get current directory to return after normalization
  rgset <- minfi::read.metharray(x, verbose = TRUE)
  mset <- minfi::preprocessIllumina(rgset,bg.correct=T)
  mset <- minfi::mapToGenome(mset)
  dataNorm<-minfi::getBeta(mset, type = "Illumina")
  return(dataNorm)
}



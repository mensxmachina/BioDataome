#' Create metadata of a GEO dataset for BioDataome
#'
#' Given a GSE id this function downloads phenotype data from GEO for a specific study, discovers control samples
#' and creates the metadata file for BioDataome
#' @param x a GSE series ID
#' @param y a GEO platform id (GPL)
#' @return a data frame of metadata with columns: sample IDs, Class and all other GEO phenotype data found in series matrix

#' @examples
#' metadata<-GSEmetadata("GSE11761","GPL570")
#' @export

GSEmetadata<-function(x,y){

  phenos<-downloadPhenotypePlatform(x,y)
  if (!is.na(phenos)){
    samples<-row.names(phenos)
    controls<-controlSamples(phenos)
    phenos<-cbind(samples=samples,class=controls[match(samples,controls[,1]),2],phenos)

  } else {
    phenos<-NA
  }

  return(phenos)

}

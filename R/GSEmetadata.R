#' Create sample phenotype metadata of a GEO dataset including control discovery
#'
#' Given a GSE id this function downloads phenotype data from GEO for a specific study, discovers control samples
#' and returns all sample phenotype metadata
#' @param x a GSE series ID
#' @param y a GEO platform id (GPL)
#' @return a data frame of metadata with columns: sample IDs, Class and all other GEO phenotype data found in series matrix

#' @examples \dontrun{
#' metadata<-GSEmetadata("GSE11761","GPL570") }
#' metadata<-BioDataome:::phenos
#' @export

GSEmetadata<-function(x,y){

  #platforms currenty curated in BioDataome
  platforms<-platformInfo$Technology

  if (missing(x))
    stop("Need to specify a GEO Series id, i.e 'GSE10026'")
  if (missing(y))
    stop("Need to specify a GEO technology i.e. 'GPL570'")
  if (!grepl("GSE[0-9]+",x))
    stop("x must be a GEO Series id, i.e 'GSE10026'")
  if (!any(grepl(y, platforms, ignore.case=TRUE)))
    stop("y must be one of the technologies: 'GPL570', 'GPL96', 'GPL6244', 'GPL1261','GPL13534'")

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

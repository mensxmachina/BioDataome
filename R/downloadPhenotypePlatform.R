#' Download series matrix from GEO for a given study and platform
#'
#' This function downloads returns the series matrix related to a given GEO Series (GSE)
#'
#' @param x a GEO Series id (GSE)
#' @param y a GEO platform id (GPL)
#' @return a data frame with the contents of the series matrix found in GEO
#' @examples \dontrun{
#' downloadPhenotypePlatform("GSE11761","GPL570") }
#' p<-BioDataome:::phenoPlatform
#' @export
#' @importFrom Biobase pData

downloadPhenotypePlatform<-function(x,y){
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


  phenos <- downloadPhenotype(x)
  if (length(phenos)>1){
    series_ind<-c()
    for (bb in 1:length(phenos)){
      series_ind[bb]<- regexpr(y,names(phenos)[bb],ignore.case = T)[1]
    }
    if (all(series_ind<0)){
      phenos<-NA
    } else {
      phenos<- Biobase::pData(phenos[[which(series_ind>0)]])
    }

  } else {
    phenos<- Biobase::pData(phenos[[1]])
  }
  phenos<-data.frame(phenos, stringsAsFactors=FALSE)
  return(phenos)
}

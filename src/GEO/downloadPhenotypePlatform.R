#' Download series matrices from GEO for a given study and platform
#'
#' This function downloads the series matrices related to a given GEO Series (GSE)
#' for the given platform
#'
#' @param x a GEO Series id (GSE)
#' @param y a GEO platform id (GPL)
#' @return a data frame with the contents of the series matrix
#' @examples
#' downloadPhenotypePlatform("GSE11761","GPL570")
#' @export

downloadPhenotypePlatform<-function(x,y){

  phenos <- downloadPhenotype(x)
  if (length(phenos)>1){
    series_ind<-c()
    for (bb in 1:length(phenos)){
      series_ind[bb]<- regexpr(y,names(phenos)[bb])[1]
    }
    phenos<- pData(phenos[[which(series_ind>0)]])
  } else {
    phenos<- pData(phenos[[1]])
  }
  return(phenos)
}

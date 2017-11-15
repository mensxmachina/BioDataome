#' Find all GSE ids for a given Entrez query
#'
#' @param x an esearch object as a result of an entrez_search query
#' @return a directory with the GSE with the compressed RAW files
#' @examples
#' downloadRaw("GSE11761")
#' @export
#' @importFrom GEOquery getGEOSuppFiles

entrezIDtoGSE<-function(x){
  count<-x$count
  etoGSE<-matrix(NA,nrow=count,ncol=2)
  for (i in 1:count){
    if (!is.na(x$ids[i])){
      a<-entrez_summary(db="gds",  id=x$ids[i])
      etoGSE[i,1]<-a$gse
      etoGSE[i,2]<-x$ids[i]
    } else {
      etoGSE[i,1]<-a$gse
      etoGSE[i,2]<-NA
    }

    print(i)
  }
  colnames(etoGSE)<-c("GSE","entrezID")
 return(etoGSE)
}


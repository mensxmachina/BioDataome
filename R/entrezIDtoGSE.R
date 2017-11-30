#' Find all GSE ids for a given Entrez query
#'
#' @param x an esearch object as a result of an entrez_search query
#' @return a matrix the first column of which is the GSE id and the second the entrezID
#' @examples
#' query GEO for all Homo sapiens studies with sample size between 200-300, measured with GPL570 and provide CEL files
#' r_search <- entrez_search(db="gds", term="Homo sapiens[ORGN] AND CEL[SFIL] AND gpl570[ACCN] AND 200:300[Number of Samples]",retmax =10000, use_history=TRUE)
#' entrezIDtoGSE(r_search)
#' @export
#' @importFrom rentrez entrez_summary

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


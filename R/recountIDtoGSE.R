#' Map a recount dataset ID to GSE ID
#'
#' @param x a recount dataset ID
#' @return a GSE series ID
#' @examples \dontrun{
#' recountIDtoGSE("SRP032775") }
#' @export
#' @importFrom rentrez entrez_search
#' @importFrom RCurl getURL


recountIDtoGSE<-function(x){
  if (missing(x))
    stop("Need to specify valid accession id, i.e 'SRP032775'")

  term<-paste0("Homo sapiens[ORGN] AND gpl11154[ACCN] AND ",x,"[All Fields]")
  r_search <- rentrez::entrez_search(db="gds", term=term,retmax =10000, use_history=TRUE)
  idd<-r_search$ids
  if(length(idd)==0){
    srpToGse<-NA
  } else {
    a<-rentrez::entrez_summary(db="gds",  id=idd)
    srpToGse<-a$accession
  }

  if( is.null(srpToGse) || is.na(srpToGse) ) {
    uri<-paste0("https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=",x)
    webpage <- RCurl::getURL(uri)
    webpage <- readLines(tc <- textConnection(webpage)); close(tc)
    pos<-grep("Center Project:",webpage)
    gse<-webpage[pos]
    if (length(gse)>=1){
      a<-unlist(strsplit(gse,"</CENTER_PROJECT_NAME></strong>"))
      srpToGse<-unlist(strsplit(a,"<CENTER_PROJECT_NAME>"))[2]
    } else {
      srpToGse<-NA
    }
  }

  if (!grepl("GSE",srpToGse)) {
    srpToGse<-NA
  }

return(srpToGse)
}


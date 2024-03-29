#' Column wise comparison of a dataset to a list of datasets.
#'
#' This function finds all datasets of the list y for which dataset x shares samples. The datasets are in the form
#' variables (probes) x samples. The number of variables (probes) in both datasets should be the same.
#' Let us assume we want to compare normalized gene expression dataset GSE86013 with datasets:
#' GSE86015, GSE9008, GSE9119. x and y can be either local paths where .Rda normalized data are stored
#' or links to the csv files in BioDataome.
#' First example runs with datasets stored in .csv in BioDataome.
#' Since these datasets are large we propose to use fread from package data.table to read datasets faster
#'
#' @param x the path to a normalized dataset x
#' @param y a character vector of all paths to datasets to compare
#' @return a character vector of all datasets for which dataset x shares at least one sample, separated by ;
#' @examples
#' x<-"http://dataome.mensxmachina.org/data/Homo%20sapiens/GPL570/GSE86013.csv"
#' y<-c("GSE86015.csv","GSE9008.csv","GSE9119.csv")
#' y<-paste0("http://dataome.mensxmachina.org/data/Homo%20sapiens/GPL570/",y)
#' commonGSEs<-compareDsetList(x,y)
#' @importFrom data.table fread
#' @export

compareDsetList<-function(x,y){

  if (length(grep("http",x))!=0) {
    d1<-data.table::fread(x, sep=",", header=T,stringsAsFactors=FALSE)
    d1<-t(d1[,2:ncol(d1)])
  } else {
    d1<-get(load(x))
  }
  t<-c()
  o <- 1
  for (i in 1:length(y)){
    if (x == y[i]){
      next
    }
    if (length(grep("http",y[i]))!=0) {
      prob.load <- try(d2<-data.table::fread(y[i], sep=",", header=T,stringsAsFactors=FALSE),TRUE)
      try(d2<-t(d2[,2:ncol(d2)]),TRUE)
    } else {
      prob.load <- try(d2<-get(load(y[i])),TRUE)
    }
    if (isTRUE(class(prob.load) == "try-error")){
      next
    }
    #--------
    if (dim(d2)[1] == dim(d1)[1]){
      t[o]<-compareDsets(d1,d2)
      o <- o + 1
    }
    #--------
  }
  #which dsets have common samples with x
  dsets<-basename(y[which(t!=0)])
  dsets<-strsplit(dsets,"[.]")
  dsets<-sapply(dsets,"[[",1)
  dsets<-paste(dsets, sep=";",collapse=";")

  return(dsets)

}

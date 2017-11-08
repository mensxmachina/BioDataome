#this function receives as input a path to a dataset x
#and a character vector of all paths to datasets to compare y
#and returns a list of GSE codes for which dataset in x has at least
#one common sample


compareDsetList<-function(x,y){
  d1<-get(load(x))
  t<-c()
  for (i in 1:length(y)){
    d2<-get(load(y[i]))
    t[i]<-compareDsets(d1,d2) 
  }
  #which dsets have common samples with x
  dsets<-basename(y[which(t!=0)])
  dsets<-unlist(strsplit(dsets,".Rda"))
  dsets<-paste(dsets, sep=";",collapse=";")
  return(dsets)
  
}

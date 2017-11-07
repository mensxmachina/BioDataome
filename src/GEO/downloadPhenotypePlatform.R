#this function receives as input a GEO Series (GSE) id, i.e."GSE11761"
#and a platform id i.e. GPL570 returns a data frame with the contents
#of the series matrix downloaded from GEO.

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
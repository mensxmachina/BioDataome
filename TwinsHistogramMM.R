options(stringsAsFactors = F)
rm(list=ls())
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/KLdiv.R')

dataRDir<-"E:/Dropbox (MXM - UoC)/Dataome/Mus Musculus/GPL1261/"
indexPath<-"E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/"
dlist<-list.files(path=dataRDir, pattern=".Rda")
pcaDir<-"E:/Dropbox (MXM - UoC)/Dataome/Mus Musculus//GPL1261-PCA/"
#list all files starting from "GSE" followed my any sequence of numbers
#and ending with .Rda
dlist<-dlist[grep(pattern="GSE[0-9]+.Rda$", dlist)]


results<-read.table(paste0(indexPath,"annotGPL1261.txt"), sep="\t", header=T, stringsAsFactors = F)
select<-which(  (results$Samples>=40) &  (results$Dublicates=="")  )

#dlist<-dlist[select]
dlist<-unlist(strsplit(dlist,".Rda"))
dlist<-dlist[which(!is.na(match(dlist,results$GSE[select])))]
dlist<-paste0(dlist,".Rda")

#dlist<-dlist[1:2]
ranksPALL<-c()
ranksQALL<-c()

for (i in 1:length(dlist)){
  d<-get(load(paste0(dataRDir,dlist[i])))
  
  ranksP<-c()
  ranksQ<-c()
  for (j in 1:10){
    if (ncol(d)==40 | ncol(d)==41 ){
      a<-20
    } else {
      a<-sample(floor(max(20,ncol(d)/4)): floor(ncol(d)/2)  ,1)  
    }
    b<-sample(1:ncol(d),a,replace=F)
    c<-setdiff(1:ncol(d),b)
    d1<-d[,b]
    d2<-d[,c]
    pca1<-prcomp(t(d1),scale=T)
    pca2<-prcomp(t(d2),scale=T)
    PQKL<-KLdiv(pca1,pca2)[4]
    rankP<-1
    rankQ<-1
    simDsetsP<-c()
    simDsetsQ<-c()
    for (k in setdiff(1:length(dlist),i)){
    pc<-get(load(paste0(pcaDir,dlist[k])))
    dd<-get(load(paste0(dataRDir,dlist[k])))
     if (  KLdiv(pca1,pc)[4] < PQKL ){
        rankQ<-min( (rankQ+1) ,10 )
        simDsetsP<-c(simDsetsP,dlist[k])
      }
      if (  KLdiv(pca2,pc)[4] < PQKL ){
        rankP<-min( (rankP+1) ,10 ) 
        simDsetsQ<-c(simDsetsQ,dlist[k])
      }
     #print(k) 
    }
    ranksP[j]<-rankP
    ranksQ[j]<-rankQ
    #print(j)  
 
  }
  ranksPALL<-cbind(ranksPALL,ranksP)
  ranksQALL<-cbind(ranksQALL,ranksQ)
print(i)
#if ( i %in% c(50,100,150,200,250,300,350,400)){
#  save(ranksPALL, file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/ranksPALLMM.Rda")
#  save(ranksQALL, file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/ranksQALLMM.Rda")
#}

}

options(stringsAsFactors = F)
rm(list=ls())

pcaDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570-PCA/"
dataRDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570/"
indexPath<-"E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/"
dlist<-list.files(path=pcaDir, pattern=".Rda")

results<-read.table(paste0(indexPath,"annotGPL570.txt"), sep="\t", header=T, stringsAsFactors = F)
combs<-t(combn(dlist,2))
explain<-list()
for (i in 1:nrow(combs)){
  d1<-get(load(paste0(pcaDir,combs[i,1])))
  d2<-get(load(paste0(pcaDir,combs[i,2])))
 
  d1<-get(load(paste0(pcaDir,"GSE4290.Rda")))
  d2<-get(load(paste0(pcaDir,"GSE13576.Rda")))
  
  data1<-t(get(load(paste0(dataRDir,"GSE4290.Rda"))))
  
  retained1<-min(which(stats:::summary.prcomp(d1)$importance[3, ]>=0.5))
  retained2<-min(which(stats:::summary.prcomp(d2)$importance[3, ]>=0.5))
  
  eigenvalues1 <- (d1$sdev)^2
  eigenvalues1 <- eigenvalues1[1:retained1]
  eigenvectors1 <-d1$rotation[,1:retained1]
  eigenvalues2 <- (d2$sdev)^2
  eigenvalues2 <- eigenvalues2[1:retained2]
  eigenvectors2 <-d2$rotation[,1:retained2]
  
  
  nvars<-nrow(eigenvectors1)
  load1<-c()
  for (i in 1:retained1){
    load1<-c(load1,order(abs(eigenvectors1[,i]))[1:50])
  }
  load1<-unique(load1)
  
  load2<-c()
  for (i in 1:retained2){
    load2<-c(load2,order(abs(eigenvectors2[,i]))[1:50])
  }
  load2<-unique(load2)
  
  all<-c(load1,load2)
  all2<-unique(all)
  
  
}
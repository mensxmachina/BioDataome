
options(stringsAsFactors = F)
rm(list=ls())
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/KLdiv.R')
pcaDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570-PCA/"
dataDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570/"

dlist<-list.files(path=pcaDir, pattern=".Rda")
d1<-get(load(paste0(pcaDir,"GSE33630.Rda"))) 
d2<-get(load(paste0(pcaDir,"GSE35570.Rda")))
obsKL<-KLdiv(d1,d2)[4]
#load dsets
dataset1<-get(load(paste0(dataDir,"GSE33630.Rda"))) 
dataset2<-get(load(paste0(dataDir,"GSE35570.Rda"))) 

KLpvalue<-function(dataset1,dataset2){
  #total number of samples
  n<-sum(dim(dataset1)[2],dim(dataset2)[2])
  #number of features
  N<-dim(dataset1)[1]
  dall<-cbind(dataset1,dataset2)
  #dall<-standardise(dall)
  n1<-dim(dataset1)[2]
  avgKL <- replicate(300, {
    b <- sample.int(n, replace=FALSE)
    newd1 <- dall[,b[1:n1] ]
    newd2 <- dall[ ,b[-c(1:n1)] ]
    d1<-prcomp(t(newd1),scale=T)
    d2<-prcomp(t(newd2),scale=T)
    return(KLdiv(d1,d2)[4])
  })
  
  #pvalue<-(sum(abs(avgKL) > abs(obsKL)) + 1) / (length(avgKL) + 1)
  return(avgKL)
}

simPath<-"E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/"

KLs<-get(load(paste0(simPath,"GPL570-405datasets-similarity-graph.Rda")))
KLs<-KLs[order(as.numeric(KLs[,10]), decreasing = F),]

tall<-get(load(paste0(simPath,"pvalues50.Rda"))) 

#tall<-c()

for (i in 51:nrow(KLs)){
 # d1<-get(load(paste0(pcaDir,KLs[i,1],".Rda"))) 
  #d2<-get(load(paste0(pcaDir,KLs[i,4],".Rda")))
  #obsKL<-KLdiv(d1,d2)[4]
  #load dsets
  dataset1<-get(load(paste0(dataDir,KLs[i,1],".Rda"))) 
  dataset2<-get(load(paste0(dataDir,KLs[i,4],".Rda"))) 
  tall<-rbind(tall,KLpvalue(dataset1,dataset2))
  print(i)
  
  
}
save(tall, file=paste0(simPath,"pvalues50.Rda"))
pvalue<-c()
for (i in 1:nrow(tall)){
#hist(tall[i,], breaks=100)
#abline(v=as.numeric(KLs[i,10]),col="red")
pvalue[i]<-(sum(abs(tall[i,]) > abs(as.numeric(KLs[i,10]))) + 1) / (ncol(tall) + 1)
}
write.table(pvalue,paste0(simPath,"pvalues50.txt"), sep="\t", row.names = F)



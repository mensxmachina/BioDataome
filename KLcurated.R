
options(stringsAsFactors = F)
rm(list=ls())
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/KLdiv.R')
pcaDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570-PCA/"
dlist<-list.files(path=pcaDir, pattern=".Rda")

#select dsets with >=40 sample size

indexPath<-"E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/"
results<-read.table(paste0(indexPath,"annotGPL570.txt"), sep="\t", header=T, stringsAsFactors = F)

#select<-which(  (results$Samples>=40) &  (results$Dublicates=="") & (results$Disease!="NA") )

#dlist<-dlist[select]

#dlist<-dlist[1:50]

combs<-t(combn(dlist,2))
KLs<-matrix(NA,nrow=nrow(combs),ncol=11)
for (i in 1:nrow(combs)){
  d1<-get(load(paste0(pcaDir,combs[i,1])))
  d2<-get(load(paste0(pcaDir,combs[i,2])))
 
  #check if they have common samples
  dsetname1<-unlist(strsplit(combs[i,1],".Rda"))
  dsetname2<-unlist(strsplit(combs[i,2],".Rda"))
  
  ##commons<-unlist(strsplit(results$Dublicates[match(dsetname1,results$GSE)],";"))
  ##iscommons<-match(dsetname2,commons)
  ##if (!is.na(iscommons)){
 ##   coms<-TRUE
 ## } else {
  ##  coms<-FALSE
  ##}
  
  #diseases
  disease1<-results$Disease[match(dsetname1,results$GSE)]
  disease2<-results$Disease[match(dsetname2,results$GSE)]
  ssize1<-results$Samples[match(dsetname1,results$GSE)]
  ssize2<-results$Samples[match(dsetname2,results$GSE)]
  
  KLs[i,]<-c(dsetname1,ssize1,disease1,dsetname2,ssize2,disease2,KLdiv(d1,d2))
  print(i)
  
}
colnames(KLs)<-c("dset1","samples1","disease1","dset2","samples2","disease2","retained1","retained2","minKL","SKL","maxKL")
write.table(KLs, file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/GPL570-949datasets-similarity-graph.txt",
            sep="\t",row.names = F, col.names = T)

save(KLs,file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/GPL570-949datasets-similarity-graph.Rda" )

#which(KLs[,6]=="TRUE")

KLs<-get(load(file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/GPL570-949datasets-similarity-graph.Rda" ))

#add if dublicates

duplicate<-c()
for (i in 1:nrow(KLs)){
  
 a<-unlist(strsplit(results$Dublicates[match(KLs[i,1],results$GSE)],";"))  
  
 duplicate[i]<-KLs[i,4]  %in%  a
 print(i) 
  
}

KLs<-cbind(KLs,duplicate)

colnames(KLs)<-c("dset1","samples1","disease1","dset2","samples2","disease2","retained1","retained2","minKL","SKL","maxKL","hasCommons")
write.table(KLs, file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/GPL570-949datasets-similarity-graph.txt",
            sep="\t",row.names = F, col.names = T)



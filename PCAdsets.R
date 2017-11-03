options(stringsAsFactors = F)
rm(list=ls())
library(Rfast)
dataRDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/RNASeq/"
indexPath<-"E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/"
writePath<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/RNASeq-PCA/"


dlist<-list.files(path=dataRDir, pattern=".Rda")
#dlist<-dlist[grep(pattern="GSE[0-9]+.Rda$", dlist)]


remain<-setdiff(dlist,list.files(writePath, pattern = ".Rda"))



for (i in remain){
  d<-get(load(paste0(dataRDir,i)))
  d<-t(d)
  #find constant columns
  sds<-colVars(d, std=F)
  cc<-which(sds<=0.0000001)
  
  noise <- runif(nrow(d),min=0,max=0.0000001) # generate the noise to add
  d[,cc] <- d[,cc] + noise  
  pcad <- prcomp(d, scale=T)
  save(pcad,file=paste0(writePath,i))
  print(i)
  
}

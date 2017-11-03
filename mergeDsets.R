options(stringsAsFactors = F)
rm(list=ls())
dataPath<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570/"
indexPath<-"E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/"
writeDir<-"E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/gene expression microarrays/Homo sapiens/GPL570/"
results<-read.table(paste0(indexPath,"annotGPL570.txt"), sep="\t", header=T)

#list all files in dataDir
dlist<-list.files(path=dataPath, pattern=".Rda")
#list all files starting from "GSE" followed my any sequence of numbers
#and ending with .Rda
dlist<-dlist[grep(pattern="GSE[0-9]+.Rda$", dlist)]
#merge all dsets

colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)


meanALL<-c()
sdALL<-c()
for (i in 1:length(dlist)){
  dset<- get(load(paste0(dataPath,dlist[i])))
  dset<-t(dset)
  #dsetsALL<-rbind(dsetsALL,dset)
  meanALL<-rbind(meanALL,colMeans(dset))
  sdALL<-rbind(sdALL,colSd(dset))
  print(i)
}
write.table(meanALL, file="E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/meansALL570.txt",
            sep="\t",row.names = F, col.names = T)
write.table(sdALL, file="E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/sdALL570.txt",
            sep="\t",row.names = F, col.names = T)

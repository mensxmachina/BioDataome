options(stringsAsFactors = F)
path<-"//algonas.csd.uoc.gr/Public/Dataome/"
setwd("E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/")
Platforms<-read.table("UsefulFiles/Platforms.txt", sep="\t", header=T)
#find gene  GEO platforms

platforms<-Platforms$Technology[grep("GPL",Platforms$Technology)]

for (i in 2:length(platforms)){
  f<-read.table(files[i],sep="\t", header=T)
  test<-paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", f$GSE)
  f2<-cbind(f,path=test)
  platform<-unique(f$Technology)
  write.table(f2, file=paste0(path,platform,".txt"),
              sep="\t",row.names = F, col.names = T)
}

#for the RNASeq
f<-read.table("//algonas.csd.uoc.gr/Public/Dataome/RNASeq.txt",sep="\t", header=T)
test<-paste0("https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=", f$GSE)
f2<-cbind(f,path=test)
platform<-unique(f$Technology)
write.table(f2, file=paste0(path,platform,".txt"),
            sep="\t",row.names = F, col.names = T)


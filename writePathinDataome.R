options(stringsAsFactors = F)
path<-"//algonas.csd.uoc.gr/Public/Dataome/"


files<-list.files(path,pattern = ".txt",full.names = T)


for (i in 1:length(files)){
  f<-read.table(files[i],sep="\t", header=T)
  test<-paste0("data/",f$Species,"/",f$Technology,"/",f$GSE,".csv")
  f2<-cbind(f,path=test)
  #f2<-f2[,1:11]
  #f2$Dublicates<-""
  platform<-unique(f$Technology)
  write.table(f2, file=paste0(path,platform,".txt"),
              sep="\t",row.names = F, col.names = T)
  print(i)
}

#to write the DOID
DOIDs<-get(load(file="E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/DOLevelsUniques.Rda"))
#DOIDs<-get(load("E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/gene expression microarrays/Dset Annotation/DOlevel0.Rda"))


for (i in 1:length(files)){
f<-read.table(files[i],sep="\t", header=T)
codes<-DOIDs[match(tolower(f$Disease),tolower(DOIDs[,1])),2]
codes2<-strsplit(codes,"DOID:")
codes2<-sapply(codes2,"[",2)
#find not NA
linksAvailable<-which(!is.na(codes2))
test<-paste0("disease-ontology.org/term/DOID%3A",codes2[linksAvailable])
toadd<-rep(NA, nrow(f))
toadd[linksAvailable]<-test
#toadd[which(is.na(codes2))]<-""
f2<-cbind(f,link=toadd)
f2<-f2[,c(1:11,13)]
f2[is.na(f2)] <- ""
platform<-unique(f$Technology)
write.table(f2, file=paste0(path,platform,".txt"),
            sep="\t",row.names = F, col.names = T)
print(i)
}

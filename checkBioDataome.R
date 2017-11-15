

dataPath<-"//algonas.csd.uoc.gr/Public/Dataome/"

Platforms<-read.table("E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/UsefulFiles/Platforms.txt", sep="\t", header=T)
i=1
indFile<-read.table(paste0(dataPath,Platforms$Technology[i],".txt"), sep="\t", header=T)

#EXISTING DATASETS IN RDA

#list all rda files

allrda<-list.files(paste0(dataPath,Platforms$Species[i],"/",Platforms$Technology[i]),".Rda")
#lsit all Annotation Rda files
annotrda<-list.files(paste0(dataPath,Platforms$Species[i],"/",Platforms$Technology[i]),"_Annot.Rda")

#find missing data rda
dataRda<-setdiff(allrda,annotrda)
gses<-unlist(strsplit(dataRda,".Rda"))
missing<-setdiff(indFile$GSE,gses)
missing

#find missing annotation rda
gses<-unlist(strsplit(annotrda,"_Annot.Rda"))
missing<-setdiff(indFile$GSE,gses)
missing

#list all csv files

allcsv<-list.files(paste0(dataPath,Platforms$Species[i],"/",Platforms$Technology[i]),".csv")
#lsit all Annotation csv files
annotcsv<-list.files(paste0(dataPath,Platforms$Species[i],"/",Platforms$Technology[i]),"_Annot.csv")

#find missing data csv
datacsv<-setdiff(allcsv,annotcsv)
gses<-unlist(strsplit(datacsv,".csv"))
missing<-setdiff(indFile$GSE,gses)
missing

#find missing annotation csv
gses<-unlist(strsplit(annotcsv,"_Annot.csv"))
missing<-setdiff(indFile$GSE,gses)
missing












#EXISTING Annotation IN RDA
annotRda<-list.files(paste0(dataPath,Platforms$Species[i],"/",Platforms$Technology[i]),"_Annot.Rda")
gses<-unlist(strsplit(annotRda,"_Annot.Rda"))
missing<-setdiff(indFile$GSE,gses)
missing
#COrrect
for (j in 172:length(missing)){
  metadata<-GSEmetadata(missing[j],Platforms$Technology[i])
  save(metadata, file=paste0(dataPath,Platforms$Species[i],"/",Platforms$Technology[i],"/",missing[j],"_Annot.Rda"))
  print(j)

}

#EXISTING data IN CSV
patt<-"GSE[0-9]+.csv$"
dataCsv<-list.files(paste0(dataPath,Platforms$Species[i],"/",Platforms$Technology[i]),patt)
gses<-unlist(strsplit(dataCsv,".csv"))
missing<-setdiff(indFile$GSE,gses)
missing

#EXISTING Annotation IN csv
annotCsv<-list.files(paste0(dataPath,Platforms$Species[i],"/",Platforms$Technology[i]),"_Annot.csv")
gses<-unlist(strsplit(annotCsv,"_Annot.csv"))
missing<-setdiff(indFile$GSE,gses)
missing






x<-paste0(dataPath,Platforms$Species[i],"/",Platforms$Technology[i],"/",gses,"_Annot.Rda")

RdaToCsv<-function(x){

d<-get(load(x))
filename<-unlist(strsplit(x,"_Annot.Rda"))

write.table(d, file=paste0(filename,"_Annot.csv"),
            sep=",",row.names = F)

}

sapply(x,RdaToCsv )
toremove<-list.files(paste0(dataPath,Platforms$Species[i],"/",Platforms$Technology[i]),"pheno", full.names = T)

file.remove(toremove)

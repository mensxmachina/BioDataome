<<<<<<< HEAD

options(stringsAsFactors = F)
library(Rfast)
library(rentrez)
library(XML)
library(RCurl)
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/diseasetoDO2.R')
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/diseasetoDO3.R')
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/compareDsets.R')
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/GSEtoDiseaseComb.R')
#this script annotates datasets
#by creating a data frame with columns:
# 1st column: GSE code
# 2nd column: Species
# 3rd column: Biological entity measured (i.e Gene expression, DNA methylation)
# 4th column: Measurement technology (i.e. GPL570)
# 5th column: Measurement type (i.e in situ oligonucleotide, high-throughput sequencing)
# 6th column: Sample size
# 7th column: Duplicates
# 8th column: disease
# 9th column: disease DO2
# 10th column: disease DO3

#Species, Biological entity measured and Measurement technology
#is given by the user

species<-"Homo sapiens"
BioEntity<-"Gene expression"
technology<-"GPL570"
type<-"in situ oligonucleotide"

#directory where data reside
dataDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570/"
#list all files in dataDir
dlist<-list.files(path=dataDir, pattern=".Rda")
#list all files starting from "GSE" followed my any sequence of numbers
#and ending with .Rda
dlist<-dlist[grep(pattern="GSE[0-9]+.Rda$", dlist)]
#remove .Rda extention
gse<-unlist(strsplit(dlist,".Rda"))


#call GSEtoDisease
#gse<-gse[1:5]
ssize<-c()
dublicates<-c()
dis<-c()
DO2<-c()
DO3<-c()
for (i in 1:length(gse)){
 dset<-get(load(paste0(dataDir,gse[i],".Rda"))) 
 ssize[i]<-dim(dset)[2]
 gse2<-setdiff(gse,gse[i])
 commons<-c()
  for (j in 1:length(gse2)){
    dset2<-get(load(paste0(dataDir,gse2[j],".Rda"))) 
    commons[j]<-compareDsets(dset,dset2)
  }
 a<-gse2[which(commons!=0)]
 dublicates[i]<-paste(a, sep=";",collapse=";")
 diseases<- GSEtoDiseaseComb(gse[i])
 if (!is.na(diseases)){
   #keep only first disease
   dis[i]<-unlist(strsplit(diseases,";"))[[1]][1]
   #Disease ontology level2
   DO2[i]<-diseasetoDO2(dis[i])
   #Disease ontology level3
   DO3[i]<-diseasetoDO3(dis[i])  
 } else {
   dis[i]<-NA
   DO2[i]<-NA
   DO3[i]<-NA
 }
 print(i) 
}
#gse2<-gse[1:3451]


#rr<-read.table(file="E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/annotGPL570NEW.txt", header=T)
results<-data.frame(GSE=gse,Species=rep(species,length(gse)),Entity=rep(BioEntity,length(gse)),
                    Technology=rep(technology,length(gse)),Type=rep(type,length(gse)),Samples=ssize,
                    Dublicates=dublicates,Disease=dis,DOlevel2=DO2,DOlevel3=DO3)
                    
#rall<-rbind(rr,results2)
write.table(results, file="E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/annotGPL570NEW.txt",
sep="\t",row.names = F, col.names = T)


results<-read.table(file="E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/annotGPL570NEW.txt")

=======
#this script annotates datasets
#by creating a data frame with columns:
# 1st column: GSE code
# 2nd column: Species
# 3rd column: Biological entity measured (i.e Gene expression, DNA methylation)
# 4th column: Measurement technology (i.e. GPL570)
# 5th column: Measurement type (i.e in situ oligonucleotide, high-throughput sequencing)
# 6th column: Sample size
# 7th column: Duplicates
# 8th column: disease
# 9th column: disease DO2
# 10th column: disease DO3

#Species, Biological entity measured and Measurement technology
#is given by the user

species<-"Homo sapiens"
BioEntity<-"Gene expression"
technology<-"GPL570"
type<-"in situ oligonucleotide"

#directory where data reside
dataDir<-"E:/Dropbox (MXM - UoC)/SCANEDSeries/Homo Sapiens/GPL570"
#list all files in dataDir
dlist<-list.files(path=dataDir, pattern=".Rda")
#list all files starting from "GSE" followed my any sequence of numbers
#and ending with .Rda
dlist<-dlist[grep(pattern="GSE[0-9]+.Rda$", dlist)]
#remove .Rda extention
gse<-unlist(strsplit(dlist,".Rda"))

#call GSEtoDisease

source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/GSEtoDiseaseComb.R')
gse2<-gse[1:20]
diseases<-sapply(gse,GSEtoDiseaseComb)
dis<-c()
for (i in 1:length(gse)){
  
 dis[i]<- GSEtoDiseaseComb(gse[i])
 print(i) 
}

write.table(dis, file="E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/disOne.txt",
sep="\t",row.names = F, col.names = F)





ssize
duplicates
>>>>>>> e877fce7ed80de3e445508d5afa39879e1d2dccb


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
technology<-"RNASeq"
type<-"expression profiling by high throughput sequencing"

#directory where data reside initially
dataDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/RNASeq/"
#list all files in dataDir
#dlist<-list.files(path=dataDir, pattern=".Rda")

#load GSE codes
ds<-read.table("G:/RNASeq/RECOUNT/srpToGse.txt")


#call GSEtoDisease
#gse<-gse[1:5]
ssize<-c()
dublicates<-c()
dis<-c()
DO2<-c()
DO3<-c()
for (i in 1:nrow(ds)){
 dset<-get(load(paste0(dataDir,ds$X[i],".Rda"))) 
 ssize[i]<-dim(dset)[2]
 gse2<-setdiff(ds$X,ds$X[i])
 commons<-c()
  for (j in 1:length(gse2)){
    dset2<-get(load(paste0(dataDir,gse2[j],".Rda"))) 
    commons[j]<-compareDsets(dset,dset2)
  }
 a<-gse2[which(commons!=0)]
 dublicates[i]<-paste(a, sep=";",collapse=";")
 diseases<- GSEtoDiseaseComb("GSE51358")
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

#if we don't have GSE try with pubmedid

length(which(!is.na(results$Disease)))


  term<-paste0("SRP010181[All Fields] AND rna seq[STRA] AND Homo sapiens[ORGN]")
r_search <- entrez_search(db="sra", term=term,retmax =10000, use_history=TRUE)
idd<-r_search$ids



results<-data.frame(GSE=ds$X,Species=rep(species,length(ds$X)),Entity=rep(BioEntity,length(ds$X)),
                    Technology=rep(technology,length(ds$X)),Type=rep(type,length(ds$X)),Samples=ssize,
                    Dublicates=dublicates,Disease=dis,DOlevel2=DO2,DOlevel3=DO3)
                    

write.table(results, file="E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/annotRNASeq.txt",
sep="\t",row.names = F, col.names = T)

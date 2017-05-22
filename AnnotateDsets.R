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
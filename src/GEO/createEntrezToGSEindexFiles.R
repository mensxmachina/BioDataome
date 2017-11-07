options(stringsAsFactors = F)
library(stringr)
library(rentrez)
source('src/GEO/entrezIDtoGSE.R', echo=TRUE)

#list of GEO platforms currently curated in BioDataome
platforms<-read.table("UsefulFiles/Platforms.txt", sep="\t", header=T)
#find gene epxression microarry GEO platforms
platforms<-platforms$Technology[which(platforms$Entity=="Gene expression" & platforms$Type=="in situ oligonucleotide")]

for (i in 1:length(platforms)){
  r_search <- entrez_search(db="gds", term=paste0("CEL[SFIL] AND ",platforms[i],"[ACCN]"),
                            retmax =10000, use_history=TRUE)
  
  eToGSE<-entrezIDtoGSE(r_search) 
  save(eToGSE, file=paste0("UsefulFiles/eToGSE_",platforms[i],".Rda"))
  
  
}



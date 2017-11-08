library(SCAN.UPC)
library(doParallel)
library(GEOquery)
library(rentrez)
options(stringsAsFactors = F)
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/src/GEO/controlSamples.R', echo=TRUE)
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/src/GEO/downloadPhenotypePlatform.R', echo=TRUE)
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/src/GEO/downloadPhenotype.R', echo=TRUE)
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/src/GEO/downloadRaw.R', echo=TRUE)
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/src/GEO/preprocessGEO.R', echo=TRUE)
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/src/GEO/compareDsetList.R', echo=TRUE)
dataPath<-"//algonas.csd.uoc.gr/Public/Dataome/"
tmpRAW<-"//algonas.csd.uoc.gr/Public/Dataome/tempRAW/"
#list of GEO platforms currently curated in BioDataome
#list of GEO platforms currently curated in BioDataome
setwd("E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/")
Platforms<-read.table("UsefulFiles/Platforms.txt", sep="\t", header=T)
#find gene epxression microarry GEO platforms
platforms<-Platforms$Technology[which(Platforms$Entity=="Gene expression" & Platforms$Type=="in situ oligonucleotide")]
species<- Platforms$Species[which(Platforms$Entity=="Gene expression" & Platforms$Type=="in situ oligonucleotide")]

#loop over all platforms
for (i in 1:length(platforms)){
  #query GEO for all series of a specific platform which also provide CEL files
  r_search <- entrez_search(db="gds", term=paste0("CEL[SFIL] AND ",platforms[i],"[ACCN]"),
                                                use_history=TRUE)
  #find number of records
  count<-r_search$count
  #call entrez_search again to get ALL entrezIDs
  r_search <- entrez_search(db="gds", term=paste0("CEL[SFIL] AND ",platforms[i],"[ACCN]"),
                            retmax =count, use_history=TRUE)
  entrezIDs<-r_search$ids
  
  #find r_search ids that are not already in BioDataome
  processedIDs<-list.files(file.path(dataPath,species[i],platforms[i]))
  processedIDs<-processedIDs[grep(pattern="GSE[0-9]+.Rda$", processedIDs)]
  processedIDs<-unlist(strsplit(processedIDs, ".Rda"))
  #load entrezID to GSE files
  eToGSE<-get(load(paste0("UsefulFiles/","eToGSE_",platforms[i],".Rda")))
  #find unique records
  eToGSEU<-paste0("GSE",unique(eToGSE[,1]))
  toProcessIDs<-setdiff(eToGSEU,processedIDs)
  
  #for all datasets not already in in BioDataome
  
  for (j in 1:length(toProcessIDs)){
    #download annotation
    phenos<-downloadPhenotypePlatform(toProcessIDs[j],platforms[i])
    samples<-row.names(phenos)
    #add control information in annotation
    controls<-controlSamples(phenos)
    phenos<-cbind(samples=samples,class=controls[match(samples,controls[,1]),2],phenos)
    #write annotation file in BioDataome
    writePath<-paste0(dataPath,species[i],"/",platforms[i],"/",toProcessIDs[j])
    write.table(phenos, file=paste0(writePath,"_Annot.csv"),
                sep=",",row.names = F)
    #download raw files
    downloadRaw(toProcessIDs[j],tmpRAW)
    untar(file.path(tmpRAW,paste0(toProcessIDs[j],"_RAW.tar")), exdir = tmpRAW)
    file.remove(file.path(tmpRAW,paste0(toProcessIDs[j],"_RAW.tar")))
    #by default downloadRaw will download all CEL files for the specific GSE
    celfiles<-list.files(tmpRAW, pattern = ".cel",ignore.case = T)
    #keep only cel files that have annotation
    celsKeep<-celfiles[startsWith(celfiles,samples)]
    #remove the rest
    toRemove<-setdiff(celfiles,celsKeep)
    file.remove(file.path(tmpRAW,toRemove))
    dataNorm<-preprocessGEO(tmpRAW,3)
    #save Rda file
    save(dataNorm, file=paste0(writePath,".Rda"))
    #save .csv file
    d<-round(t(dataNorm), digits = 4)
    #create proper column for csv
    row.names(d)<-samples
    d<-cbind.data.frame(samples=samples,d)
    # write data in csv
    write.table(d, file=paste0(writePath,".csv"), 
                sep=",",row.names = F)
    #clear tmpRAW
    file.remove(dir(tmpRAW, full.names=TRUE)) 
    #update index file
    
    #find any datasets that share samples 
    x<-paste0(writePath,".Rda")
    y<-paste0(dataPath,species[i],"/",platforms[i],"/",processedIDs,".Rda")
    commonGSEs<-compareDsetList(x,y)
    
    #Annotate dset wiht Disease Ontology terms
    
    
    
    
    
    }
}
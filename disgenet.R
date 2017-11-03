options(stringsAsFactors = F)
rm(list=ls())
library(Biobase)
library(GEOquery)

#Download GPL file, put it in the current directory, and load it:
gpl570 <- getGEO('GPL570', destdir=".")
geneSymbols<-Table(gpl570)[,c("ID","Gene Symbol")]

simPath<-"E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/"
pcaDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570-PCA/"
profile<-get(load(paste0(simPath,"GPL570-949datasets-profile.Rda")))

select<-which(profile[,3] =="chronic obstructive pulmonary disease" )

#find involved
dsets<-profile[select,1]
genesALL<-c()
for (i in 1:length(dsets)){
  d<-get(load(paste0(pcaDir,dsets[i],".Rda")))
  plot(d$x[,1],d$x[,2])
  
  
  
  
  cP<-min(which(stats:::summary.prcomp(d)$importance[3, ]>=0.5))
  probe<-c()
  for (j in 1:cP){
    Pmax<-d$rotation[,j] 
    involved1.5<-which( abs(Pmax) >4/(sqrt(length(Pmax))) )
    #involved1.5<-which( abs(Pmax) >5*(mean(abs(Pmax))) )
    #involved1.5<-which( abs(Pmax) >quantile(abs(Pmax),0.99) )
    #take probe name of max loading to this PC
    probe<-c(probe, names(Pmax)[involved1.5])
  }
  
  
  #Pmax<-d$rotation[,1]
  #a<-sort(abs(Pmax), decreasing = T)
  #involved1.5<-which( abs(Pmax) >1.5/(sqrt(length(Pmax))) )
  #involved1.5<-which( abs(Pmax) >quantile(abs(Pmax),0.99) )
  #involved1.5<-which( abs(Pmax) >2*(mean(a)) )
  #involved1.5<-a[1:10]
  #probes<-names(Pmax)[involved1.5]
  #probes<-names(Pmax)[1:10]
  
  genesALL<-c(genesALL,  geneSymbols$`Gene Symbol`[which(geneSymbols$ID %in% probe)]  )
  print(i)
  
}

genes<-unique(genesALL)
#geneSymbolsALL<-unlist(strsplit(geneSymbols$`Gene Symbol`," /// "))
#geneSymbolsALL<-unique(geneSymbolsALL)

genes<-unlist(strsplit(genes," /// "))
genesAIS<-read.table(file="E:/Dropbox (MXM - UoC)/Documents/Database paper/genesCOPD.txt")

length(which(genes %in% genesAIS$V1))
length(intersect(genes,genesAIS$V1))


length(which(genes  %in% geneSymbolsALL))


length(intersect(genes,geneSymbolsALL))
genes[which(genes %in% genesAIS$V1)]

which(genes %in% "CTLA4")
which(!is.na((match("IFNA1",genes))))

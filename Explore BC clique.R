options(stringsAsFactors = F)
rm(list=ls())
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/explainKL.R')
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/explainKLonly.R', echo=TRUE)

pcaDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570-PCA/"

BCClique<-c("GSE20685","GSE26639","GSE54002","GSE43365", "GSE71258",
            "GSE58812", "GSE42568", "GSE58984", "GSE16446", "GSE43358","GSE48906")

kappa<-c(10,50,100,500, 1000)

k2<-seq(2000,54000,2000)
kappa<-c(kappa,k2)

geneListALL<-list()
for (j in 1:length(kappa)){
  geneList<-list()
  
  for (i in 1:(length(BCClique)-2)){
    dsets<-BCClique[1:(i+2)]
    dsets<-paste0(dsets,".Rda")
    geneList[[i]]<-explainKLonly(dsets,kappa[j])  
    
    #print(i)
    
  }
  
  geneListALL[[j]]<-geneList
  
  
 print(j) 
  
}

aa<-c()
for (j in 1:length(kappa)){
geneList<-geneListALL[[j]]
jaccALL<-c()
for (i in 1:(length(geneList)-1)){
  jaccALL[i]<-length(intersect(geneList[[i]],geneList[[i+1]])) / length( union(geneList[[i]],geneList[[i+1]]) )
  
}

aa<-rbind(aa,jaccALL)
print(j)
}

dsets4<-c("GSE20685","GSE26639","GSE54002","GSE43365")
dsets4<-paste0(dsets4,".Rda")
dsets<-paste0(BCClique,".Rda")
tt<-explainKLonly(dsets,6000)  

Symbol <- getSYMBOL(row.names(d1$rotation), "hgu133plus2.db")
genes<-Symbol[tt]
genes2<-genes[!is.na(genes)]
genes2<-unique(genes2)
write.table(genes2,file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/explainKL/Breast Cancer/cliqueBC510genes.txt",
            sep="\t", row.names = F, quote = F, col.names = F)






source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/explainKLonly.R', echo=TRUE)
library("annotate")
library('hgu133plus2.db')
d1<-get(load(paste0(pcaDir,"GSE26639.Rda")))

entrez<-getEG(row.names(d1$rotation), "hgu133plus2.db")
library(ReactomePA)


dsets<-c("GSE20685","GSE26639","GSE54002","GSE43365", "GSE71258",
         "GSE58812", "GSE42568", "GSE58984", "GSE16446", "GSE43358","GSE48906")

dsets<-c("GSE20685","GSE26639","GSE54002","GSE43365","GSE71258","GSE58812","GSE42568")


dsets<-paste0(dsets,".Rda")
kappa<-c(10,50,100,500, 1000)

k2<-seq(2000,54000,2000)
kappa<-c(kappa,k2)
paths<-c()
for (i in kappa){
  #dsets<-dsetsC[[6]][i,]
  tt<-explainKLonly(dsets,i)
  
  tt<-unlist(tt)
  
  entrezGenes<-entrez[tt]
  entrezGenes2<-entrezGenes[!is.na(entrezGenes)]
  entrezGenes2<-unique(entrezGenes2)
  kk1 <- enrichPathway(entrezGenes2, organism = "human",pvalueCutoff=0.05, readable=T)
  if (length(kk1)!=0){
    #length(entrezGenes2)
    paths<-c(paths,length(kk1@result$ID))
  } else {
    paths<-c(paths,0)
  }
  
  
  
}
#pathsALL<-c()
pathsALL<-rbind(pathsALL,paths)



kappa<-c(10,50,100,500, 1000)

k2<-seq(2000,54000,2000)
kappa<-c(kappa,k2)
jacs<-c()
for (i in kappa){
  #dsets<-dsetsC[[6]][i,]
  jacs<-c(jacs,explainKL(dsets,i)[2])
 
  print(i)
  
}

#jacsall<-c()
jacsall<-rbind(jacsall,jacs)









Symbol <- getSYMBOL(row.names(d1$rotation), "hgu133plus2.db")
genes<-Symbol[intersections[[length(intersections)]]]
genes2<-genes[!is.na(genes)]
genes2<-unique(genes2)
write.table(genes2,file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/explainKL/Breast Cancer/cliqueBC11dsets6000.txt",
            sep="\t", row.names = F, quote = F, col.names = F)



entrezGenes<-entrez[intersections[[length(intersections)]]]
entrezGenes2<-entrezGenes[!is.na(entrezGenes)]
entrezGenes2<-unique(entrezGenes2)


kk1 <- enrichPathway(entrezGenes2, organism = "human",pvalueCutoff=0.05, readable=T)




kk1@result$Description
kk1@result$Description




library("pathview")
#thyroid
#05216

#breast
#05224

hsa05224 <- pathview(gene.data  = entrezGenes2,
                     pathway.id = "hsa05224",
                     species    = "hsa")

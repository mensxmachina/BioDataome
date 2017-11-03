options(stringsAsFactors = F)
rm(list=ls())
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/explainKL.R')
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/cliqueKL.R')
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/KLdiv.R')
source('E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/DBIR/KLparams.R', echo=TRUE)
pcaDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570-PCA/"
dataDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570/"
indexPath<-"E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/"
dlist<-list.files(path=pcaDir, pattern=".Rda")

results<-read.table(paste0(indexPath,"annotGPL570.txt"), sep="\t", header=T, stringsAsFactors = F)
load("E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/GPL570-405datasets-similarity-graph.Rda")
dsets<-read.table(paste0("E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/405dsets.txt"), sep="\t", header=T)

d1<-get(load(paste0(pcaDir,"GSE26639.Rda")))


KLs<-as.data.frame(KLs)
KLs[,7]<-as.numeric(KLs[,7])
KLs[,8]<-as.numeric(KLs[,8])

KLsNew<-KLs[ which(KLs[,7]>=3),]

KLsNew<-KLsNew[ which(KLsNew[,8]>=3),]
library(pracma)
df <- as.matrix(squareform(as.numeric(KLsNew$SKL)) )
row.names(df)<-unique(c(KLsNew$dset1,KLsNew$dset2))
colnames(df)<-row.names(df)


#breast cancer low KL
dsetsC<-c("GSE26639","GSE54002","GSE20685")
#dsetsC<-c("GSE26639.Rda","GSE54002.Rda","GSE20685.Rda")
#thyroid cancer low KL
dsetsC<-c("GSE35570","GSE33630","GSE60542","GSE29265")


dsetsC<-c("GSE35571","GSE18123","GSE33943",
          "GSE9960","GSE69961","GSE39088",
          "GSE61629","GSE57065","GSE6575",
          "GSE37171")

dsetsC<-c("GSE13294","GSE37364","GSE26682",
          "GSE9348","GSE64857","GSE28702")

dsetsC<-c("GSE13294","GSE37364","GSE26682",
          "GSE9348")
#lung cancer 1
dsetsC<-c("GSE30219","GSE43580","GSE50081")
#lung cancer 2
dsetsC<-c("GSE40791","GSE31210","GSE27262")

#lung cancer 3-clique
dsetsC<-c("GSE43580","GSE31210","GSE30219","GSE50081")
#lung big
dsetsC<-c("GSE30219","GSE40791","GSE31210","GSE50081","GSE27262","GSE10245","GSE18842","GSE43580")


dsetsC<-c("GSE12992","GSE37418","GSE10327")
avC<-cliqueKL(dsetsC,df)
dsetsC<-paste0(dsetsC,".Rda")

#randomly select datasets the same size as the clique
explainClique<-c()
explainRandomGJ<-c()
explainRandom<-c()
explainCliqueGJ<-c()
avCR<-c()
kappa<-c(2,5,7,10,50,100,500, 1000)

k2<-seq(1000,54000,2000)
kappa<-c(kappa,k2)
#kappa<-c(2,5,7)
#kappa<-seq(38000,54000,2000)
for (i in kappa){
  dsets<-sample(setdiff(unique(c(KLsNew$dset1,KLsNew$dset2)),dsetsC),length(dsetsC),replace = F)
  avCR<-c(avCR,cliqueKL(dsets,df))
  dsets<-paste0(dsets,".Rda")
  explainRandom<-c(explainRandom,explainKL(dsets,i)[1] )
  explainClique<-c(explainClique,explainKL(dsetsC,i)[1])
  
  explainRandomGJ<-c(explainRandomGJ,explainKL(dsets,i)[2] )
  explainCliqueGJ<-c(explainCliqueGJ,explainKL(dsetsC,i)[2])
  
  
  print(i)
}

#plot results

dat<-cbind(explainRandom,explainClique,explainRandomGJ,
           explainCliqueGJ)
matplot(dat, type = c("b"),pch=1,col = 1:4,
        xaxt="n", main ="Thyroid Cancer-4dsets",
        xlab="number of variables", ylab="Jaccard") #plot
axis(1, at = c(1:length(kappa)), labels=kappa , las=1)
legend("topleft", legend = c("Random","Clique", "RandomGJ", "CliqueGJ"), col=1:4, pch=1) # optional legend


#study k for lung cancers

#lung big
dsetsC<-c("GSE30219","GSE40791","GSE31210","GSE50081","GSE27262","GSE10245","GSE18842","GSE43580",
          "GSE19188","GSE10445", "GSE12667", "GSE33532", "GSE77803")
dsetsC<-paste0(dsetsC,".Rda")

combs<-t(combn(dsetsC,3))
kappa<-c(2,5,7,10,50,100,500, 1000)
k2<-seq(1000,20000,2000)
kappa<-c(kappa,k2)

#explainClique<-c()
explainRandomGJ<-c()
#explainRandom<-c()

findK<-c()
for (j in 47:nrow(combs)){
  dsetsC<-combs[j,]
  #kappa<-seq(38000,54000,2000)
  explainCliqueGJ<-c()
  for (i in kappa){
    #dsets<-sample(setdiff(unique(c(KLsNew$dset1,KLsNew$dset2)),dsetsC),length(dsetsC),replace = F)
    #avCR<-c(avCR,cliqueKL(dsets,df))
    #dsets<-paste0(dsets,".Rda")
    #explainRandom<-c(explainRandom,explainKL(dsets,i)[1] )
    #explainClique<-c(explainClique,explainKL(dsetsC,i)[1])
    
    #explainRandomGJ<-c(explainRandomGJ,explainKL(dsets,i)[2] )
    explainCliqueGJ<-c(explainCliqueGJ,explainKL(dsetsC,i)[2])
    
    
    #print(i)
  } 
  
  findK<-rbind(findK,explainCliqueGJ)
 print(j) 
  
}


matplot(t(findK), type = c("b"),pch=1,
        xaxt="n", main ="Lung Cancer-3dsets-combinations",
        xlab="number of variables", ylab="Jaccard") #plot
axis(1, at = c(1:length(kappa)), labels=kappa , las=1)
#legend("topleft", legend = c("Random","Clique", "RandomGJ", "CliqueGJ"), col=1:4, pch=1) # optional legend


#translate to genes

library("annotate")
library('hgu133plus2.db')
entrez<-getEG(row.names(d1$rotation), "hgu133plus2.db")
#entrez<-entrez[!is.na(entrez)]11
Symbol <- getSYMBOL(row.names(d1$rotation), "hgu133plus2.db")
#Symbol<-Symbol[!is.na(Symbol)]
library(ReactomePA)

genes<-Symbol[intersections[[length(intersections)]]]
entrezGenes<-entrez[intersections[[length(intersections)]]]
kk1 <- enrichPathway(entrezGenes[!is.na(entrezGenes)], organism = "human",pvalueCutoff=0.05, readable=T)
#randmo
genesH<-genes[!is.na(genes)]

paths<-c()
for (i in 1:1000){
  entrezGenesRandom<-sample(entrez,length(intersections[[2]]), replace = F)
  kk2 <- enrichPathway(entrezGenesRandom[!is.na(entrezGenesRandom)], organism = "human",pvalueCutoff=0.05, readable=T)
  paths[i]<-length(kk2@result$ID)
 print(i) 
}

length(kk1@result$ID)

head(summary(kk1))
barplot(kk1, showCategory=73)

dotplot(kk1, showCategory=15)

library("pathview")
#thyroid
#05216

#breast
#05224

hsa05216 <- pathview(gene.data  = entrezGenes[!is.na(entrezGenes)],
                     pathway.id = "hsa05216",
                     species    = "hsa")
library(DOSE)
x <- enrichDO(gene          = kk1,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = entrez,
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE) 
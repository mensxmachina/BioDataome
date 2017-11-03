library(igraph)
options(stringsAsFactors = F)
rm(list=ls())

load("E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/GPL570-405datasets-similarity-graph.Rda")
indexPath<-"E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/"
results<-read.table(paste0(indexPath,"annotGPL570.txt"), sep="\t", header=T)
dsets<-read.table(paste0("E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/405dsets.txt"), sep="\t", header=T)

pcaDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570-PCA/"

#nodes2<-read.csv("C:/Users/kliolak/Downloads/netscix2016/Dataset1-Media-Example-NODES.csv", header=T, as.is=T)

#links2 <-read.csv("C:/Users/kliolak/Downloads/netscix2016/Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)
#select those with more than 2 retained eigenvalues
KLs<-as.data.frame(KLs)
KLs[,7]<-as.numeric(KLs[,7])
KLs[,8]<-as.numeric(KLs[,8])

KLsNew<-KLs[ which(KLs[,7]>=3),]


KLsNew<-KLsNew[ which(KLsNew[,8]>=3),]

links<-data.frame(from=KLsNew[,1], to=KLsNew[,4], weight=KLsNew[,10])
links$weight <- as.numeric(links$weight)
links$from <- as.character(links$from)
links$to <- as.character(links$to)



linksSorted<-links[order(links$weight),]
av<-mean(links$weight)
bv<-which(links$weight<=av)
links<-linksSorted[1:av,]
#links<-linksSorted
nodesUsed<-unique(c(links$from,links$to))
nodes<-results[match(nodesUsed,results$GSE), c(1,8,9)]

names(nodes)<-c("id","disease", "diseaseLevel")  
#nodes<-cbind(nodes,group1=nodes$disease, group2=nodes$diseaseLevel)
#library(visNetwork)
#visNetwork(nodes, links, height = "700px", width = "100%") %>%
#  visOptions(selectedBy = "group", 
#             highlightNearest = TRUE, 
#             nodesIdSelection = TRUE) %>%
#  visPhysics(stabilization = T,repulsion=list(springLength=1000,springConstant=0)) %>%
  #  visEdges(length=400) %>%
#  visLegend(width=0.2) 



net <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
E(net)$weight=as.numeric(links$weight)

#write_graph(net, "E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/forCyto4000.gml",
#           format = "gml") 
k=20
dsetsC<-list()
test<-cliques(net, min = k)
while ( (length(cliques(net, min = k))!=0) & (k>2)    ){
  test<-cliques(net, min = k) 
  dsetsC[[k]]<-t(matrix(names(unlist(test)), nrow=k, ncol=length(test) ))
  k=k-1
}

subset<-c()
for (i in 1:nrow(dsetsC[[4]])){
  for (j in 1:nrow(dsetsC[[5]])){
    
    if (length(dsetsC[[4]][i,] %in% dsetsC[[5]][j,])==length(dsetsC[[4]][i,])) {
      subset<-c(subset,i)
    }
  }
  
  
}







library("annotate")
library('hgu133plus2.db')
d1<-get(load(paste0(pcaDir,"GSE26639.Rda")))

entrez<-getEG(row.names(d1$rotation), "hgu133plus2.db")
library(ReactomePA)




kappa<-c(10,50,100,500, 1000)

k2<-seq(1000,54000,2000)
kappa<-c(kappa,k2)
paths<-matrix(NA,nrow=nrow(dsetsC[[6]]), ncol=length(kappa) )
pathsR<-matrix(NA,nrow=nrow(dsetsC[[6]]), ncol=length(kappa) )
for (i in 1:nrow(dsetsC[[6]])){
  dsets<-dsetsC[[6]][i,]
  dsets<-paste0(dsets,".Rda")
  pathsL<-c()
  pathsLR<-c()
  for (j in kappa){
    tt<-explainKLonly(dsets,j)
    
    tt<-unlist(tt)
    
    entrezGenes<-entrez[tt]
    entrezGenes2<-entrezGenes[!is.na(entrezGenes)]
    entrezGenes2<-unique(entrezGenes2)
    kk1 <- enrichPathway(entrezGenes2, organism = "human",pvalueCutoff=0.05, readable=T)
    if (length(kk1)!=0){
      #length(entrezGenes2)
      pathsL<-c(pathsL,length(kk1@result$ID))
    } else {
      pathsL<-c(pathsL,0)
    }
      #take random genes
      entrezGenesR<-entrez[sample(1:length(entrez),length(tt))]
      entrezGenes2R<-entrezGenesR[!is.na(entrezGenesR)]
      entrezGenes2R<-unique(entrezGenes2R)
      entrezGenes2R<-entrezGenes2R[1:length(entrezGenes2)]
      kk2 <- enrichPathway(entrezGenes2R, organism = "human",pvalueCutoff=0.05, readable=T)
      if (length(kk2)!=0){
        pathsLR<-c(pathsLR,length(kk2@result$ID))   
      } else {
        pathsLR<-c(pathsLR,0) 
        
      }
      
      
      #length(entrezGenes2R)
      
      
    
  }
  paths[i,]<-pathsL
  pathsR[i,]<-pathsLR
  print(i)  
}


matplot(t(paths), type = c("b"),pch=1,
        xaxt="n",
        xlab="Genes", ylab="Pathways") #plot
axis(1, at = c(1:length(kappa)), labels=kappa , las=1)

matplot(t(pathsR), type = c("b"),pch=1,
        xaxt="n",
        xlab="Genes", ylab="Pathways") #plot
axis(1, at = c(1:length(kappa)), labels=kappa , las=1)






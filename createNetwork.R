options(stringsAsFactors = F)
rm(list=ls())
library(visNetwork)
library(reshape)

simPath<-"E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/"

KLs<-get(load(paste0(simPath,"GPL13534-117datasets-similarity-graph.Rda")))

#select rows with eigenvalues >1

select<-which(  (KLs[,7]>=4) & (KLs[,8]>=4) )
KLs<-KLs[select,]



myedges <- as.data.frame(KLs[,c(1,4,10)],stringsAsFactors = F)
names(myedges) <- c("from", "to", "weight")
myedges.sorted<-myedges[order(myedges[,3], decreasing = F),]

indexPath<-"E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/"
results<-read.table(paste0(indexPath,"annotGPL13534.txt"), sep="\t", header=T, stringsAsFactors = F)

#select<-which(  (results$Samples>=40) &  (results$Dublicates=="") & (results$Disease!="NA") )

select<-which(  (results$Samples>=40) &  (results$Dublicates=="") )

nodesALL<-results[select,c("GSE","Disease")]
mynodes<-cbind(id=nodesALL$GSE, label=nodesALL$Disease ,group=nodesALL$Disease) 
mynodes<-as.data.frame(mynodes,stringsAsFactors=F )
nnumber<-20
#links.pruned<-links.sorted[1:nnumber,]

myedges.pruned<-myedges.sorted[1:nnumber,]

#update mynodes
nodes.left<-unique(c(myedges.pruned[,1],myedges.pruned[,2]))
#mynodes.sorted<-arrange(mynodes, organID)
mynodes.pruned<-mynodes[match(nodes.left,mynodes[,1]),]

#mynodes.pruned<-cbind(mynodes[match(nodes.left,mynodes[,1]),],shape="dot")
visNetwork(mynodes.pruned, myedges.pruned, height = "700px", width = "100%") %>%
  visOptions(selectedBy = "group", 
             highlightNearest = TRUE, 
             nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = T,repulsion=list(springLength=1000,springConstant=0)) %>%
#  visEdges(length=400) %>%
  visLegend(width=0.2) 



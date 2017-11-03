options(stringsAsFactors = F)
rm(list=ls())
ranksPALL<-get(load(file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/ranksPALL.Rda"))
ranksQALL<-get(load(file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/ranksQALL.Rda"))

write.table(ranksPALL, file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/GPL1261-twins-histogram-P.txt",
            sep="\t",row.names = F, col.names = T)
write.table(ranksQALL, file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/GPL1261-twins-histogram-Q.txt",
            sep="\t",row.names = F, col.names = T)

write.table(dlist, file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/405dsets.txt",
            sep="\t",row.names = F, col.names = T)

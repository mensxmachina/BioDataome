#see ggeo note 16, 1st task

options(stringsAsFactors = F)
rm(list=ls())

dataDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570/"
dset<-get(load(paste0(dataDir,"GSE31162.Rda")))
#pca <- prcomp(t(dset), scale=T)
#eigenvalues<-(pca$sdev)^2

step=20
steps<-floor(ncol(dset)/step)
rest<-c(1:ncol(dset))
ds<-c()
results<-matrix(NA, nrow = steps, ncol=8)
for (i in 1:steps){
  #choose 20
  d<-sample(rest,step, replace = F)
  
  ds<-cbind(ds,dset[,d])
  #keep the rest
  rest<-setdiff(rest,d)
  
  pcad <- prcomp(t(ds), scale=T)
  eigens<-(pcad$sdev)^2
  retained<-min(which(stats:::summary.prcomp(pcad)$importance[3, ]>=0.5))
  #noise turnpoint
  ind<-c()
  for (k in 1:(length(eigens)-1)){
    if ( sd(eigens[k:length(eigens)])<= mean(eigens[k:length(eigens)]) ){
      ind<-k  
      break 
    }
  }
  lnoise<-mean(eigens[ind:length(eigens)])  
  sdnoise<-sd(eigens[ind:length(eigens)])
  snr<-eigens[1]/lnoise
  
  results[i,]<-c( eigens[1],eigens[2],eigens[3], retained, ind, lnoise,sdnoise, snr )
  print(i)

}

colnames(results)<-c("l1","l2","l3","retained","noiseTurn","lnoise","sdnoise","snr")
write.table(results, file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/GPL570-GSE31162-stepwise-spectrum.txt",
            sep="\t",row.names = F, col.names = T)

save(profile,file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/GPL570-GSE31162-stepwise-spectrum.Rda" )


#entire dset

pcad <- prcomp(t(dset), scale=T)
eigens<-(pcad$sdev)^2
retained<-min(which(stats:::summary.prcomp(pcad)$importance[3, ]>=0.5))
#noise turnpoint
ind<-c()
for (k in 1:(length(eigens)-1)){
  if ( sd(eigens[k:length(eigens)])<= mean(eigens[k:length(eigens)]) ){
    ind<-k  
    break 
  }
}
lnoise<-mean(eigens[ind:length(eigens)])  
sdnoise<-sd(eigens[ind:length(eigens)])
snr<-eigens[1]/lnoise




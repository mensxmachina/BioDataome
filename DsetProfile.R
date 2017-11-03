#see ggeo note 16, 2nd task

options(stringsAsFactors = F)
rm(list=ls())

pcaDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570-PCA/"
dataRDir<-"E:/Dropbox (MXM - UoC)/Dataome/Homo Sapiens/GPL570/"
indexPath<-"E:/Dropbox (MXM - UoC)/Documents/Data-Based IR/Paper Experiments/Annotation/"
dlist<-list.files(path=pcaDir, pattern=".Rda")

results<-read.table(paste0(indexPath,"annotGPL570.txt"), sep="\t", header=T, stringsAsFactors = F)

profile<-matrix(NA,nrow=length(dlist),ncol=23)
for (i in 1:length(dlist)){
  pca<-get(load(paste0(pcaDir,dlist[i])))
  eigenvalues <- (pca$sdev)^2
  eigenvectors <-pca$rotation
  retained<-min(which(stats:::summary.prcomp(pca)$importance[3, ]>=0.5))
  #noise turnpoint
  ind<-c()
  for (k in 1:(length(eigenvalues)-1)){
    if ( sd(eigenvalues[k:length(eigenvalues)])<= mean(eigenvalues[k:length(eigenvalues)]) ){
      ind<-k  
      break 
    }
  }
  
  if (!is.null(ind)){
    lnoise<-mean(eigenvalues[ind:length(eigenvalues)])  
    sdnoise<-sd(eigenvalues[ind:length(eigenvalues)]) 
  } else {
    ind<-NA
    lnoise<-NA  
    sdnoise<-NA
  }
  
  snr<-eigenvalues[1]/lnoise
  Pmax<-eigenvectors[,1]
  minAngle<-min(acos(abs(Pmax)))*180/pi
  involved1.5<-length( which( abs(Pmax) >1.5/(sqrt(length(Pmax))) ))
  involved2.0<-length( which( abs(Pmax) >2.0/(sqrt(length(Pmax))) ))
  involved2.5<-length( which( abs(Pmax) >2.5/(sqrt(length(Pmax))) ))
  involved3.0<-length( which( abs(Pmax) >3.0/(sqrt(length(Pmax))) ))
  involved3.5<-length( which( abs(Pmax) >3.5/(sqrt(length(Pmax))) ))
  involved4.0<-length( which( abs(Pmax) >4.0/(sqrt(length(Pmax))) ))
  
  snrEstim<-snr+(nrow(eigenvectors)-ncol(eigenvectors))*0.25
  delta<-atan( 1/( sqrt(snrEstim)*sqrt(ncol(eigenvectors)) ) )*180/pi
  theta<-asin( (1/sqrt(snrEstim)) * sqrt(nrow(eigenvectors)/ncol(eigenvectors))  )*180/pi
  
  profile[i,]<-c( results$GSE[i], ncol(eigenvectors),results$Disease[i], eigenvalues[1:5],
                  retained, ind,lnoise,sdnoise,minAngle,involved1.5,involved2.0,involved2.5,
                  involved3.0,involved3.5,involved4.0,snrEstim,snr,delta,theta)

  print(i)  
}


colnames(profile)<-c("GSE","Samples", "disease","1st eigenvalue", "2nd eigenvalue","3rd eigenvalue",
                     "4th eigenvalue","5th eigenvalue","retained","noiseTurn","lnoise","sdnoise",
                     "minAngle", "involved-1.5","involved-2.0","involved-2.5","involved-3.0","involved-3.5",
                     "involved-4.0","SNR-estimation","snr","delta","theta")
write.table(profile, file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/GPL570-949datasets-profile.txt",
            sep="\t",row.names = F, col.names = T)

save(profile,file="E:/Dropbox (MXM - UoC)/Paper - Landscape of Human Dataome - 2017/Results/GPL570-949datasets-profile.Rda" )

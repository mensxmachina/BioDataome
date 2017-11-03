
explainKLonly2dsets<-function(dsets,k){
#dsets is a character vector of three dsets
#i.e dsetsC<-c("GSE26639.Rda","GSE54002.Rda","GSE20685.Rda")
#k is the number of genes we want to select


  d1<-get(load(paste0(pcaDir,dsets[1])))
  d2<-get(load(paste0(pcaDir,dsets[2])))
  
  # d1<-get(load(paste0(pcaDir,"GSE4290.Rda")))
  #d2<-get(load(paste0(pcaDir,"GSE13576.Rda")))
  
  
  
  retained1<-min(which(stats:::summary.prcomp(d1)$importance[3, ]>=0.5))
  retained2<-min(which(stats:::summary.prcomp(d2)$importance[3, ]>=0.5))
  
  eigenvalues1 <- (d1$sdev)^2
  eigenvalues1 <- eigenvalues1[1:retained1]
  eigenvectors1 <-d1$rotation[,1:retained1]
  eigenvalues2 <- (d2$sdev)^2
  eigenvalues2 <- eigenvalues2[1:retained2]
  eigenvectors2 <-d2$rotation[,1:retained2]
  
  
  nvars<-nrow(eigenvectors1)
  
  
  alpha<-1/(nvars-retained1) * (nvars - sum(eigenvalues1[1:retained1]) )
  
  beta<-1/(nvars-retained2) * (nvars - sum(eigenvalues2[1:retained2]) )
  
  coefs<-  ( eigenvalues1/alpha -1  ) %*% ( 1- beta/t(eigenvalues2 )) + 
    (1- alpha/eigenvalues1) %*% ( t(eigenvalues2)/beta - 1)
  
  
  s = rep(1, nvars)
  t = rep(0, nvars)
  
  iter <- 0
  maxiters <- 150

  # B<-c()
  for (j in 1:nvars){
    
    if (( s[j]!=t[j]   & iter<=maxiters   )){
      iter<-iter+1
      print( paste0('Iteration ', iter )    )
      
      t=s
      
      W = coefs * (t(eigenvectors1[t,]) %*% eigenvectors2[t,])
      rowsums = rowSums( (eigenvectors1 %*% W ) * eigenvectors2  )     
      
      indexes<-order(rowsums, decreasing = T)
      
      s = logical(nvars)
      s[indexes[1:k]] = TRUE
    stopifnot(iter != maxiters)
      
      
    }
    
  }
  
  explainKL<- which(s)


return(explainKL)
}

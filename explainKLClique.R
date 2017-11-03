
explainKLClique<-function(dsets,k){
#dsets is a character vector of three dsets
#i.e dsetsC<-c("GSE26639.Rda","GSE54002.Rda","GSE20685.Rda")
#k is the number of genes we want to select


combs1<-t(combn(dsets,2))
KLsALL<-c()
KLsALLRandom<-c()

s = rep(1, nvars)
t = rep(0, nvars)

iter <- 0
maxiters <- 150
# B<-c()
#KL<-list()
for (j in 1:nvars){
  
  if (( s[j]!=t[j]   & iter<=maxiters   )){
    iter<-iter+1
    print( paste0('Iteration ', iter )    )
    
    t=s
    
    rowsumsALL<-c()
    KLs<-c()
    KLsRandom<-c()
    for (i in 1:nrow(combs1)){
      d1<-get(load(paste0(pcaDir,combs1[i,1])))
      d2<-get(load(paste0(pcaDir,combs1[i,2])))
      
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
    
    
      W = coefs * (t(eigenvectors1[t,]) %*% eigenvectors2[t,])
      rowsums = rowSums( (eigenvectors1 %*% W ) * eigenvectors2  )     
      rowsumsALL<-cbind(rowsumsALL,rowsums)
      
      m<-which(t==1)
      mRandom<-sample(1:nvars, k, replace=F)
      d1<-get(load(paste0(dataDir,combs1[i,1])))
      d11<-d1[m,]
      d11<-t(d11)
      pc1<-prcomp(d11, scale=T)
      
      d3<-d1[mRandom,]
      d3<-t(d3)
      pc3<-prcomp(d3, scale=T)
      
      
      d2<-get(load(paste0(dataDir,combs1[i,2])))
      d22<-d2[m,]
      d22<-t(d22)
      pc2<-prcomp(d22, scale=T)
      
      d4<-d2[mRandom,]
      d4<-t(d4)
      pc4<-prcomp(d4, scale=T)
      
      KLs<-c(KLs,KLdiv(pc1,pc2)[4])
      KLsRandom<-c(KLsRandom,KLdiv(pc3,pc4)[4])
      }
    rowsums<-rowSums(rowsumsALL)
    indexes<-order(rowsums, decreasing = T)
    
    s = logical(nvars)
    s[indexes[1:k]] = TRUE
    
    #KLsALL<-rbind(KLsALL,KLs)
    
    stopifnot(iter != maxiters)

  }
 
  
}
#KL[[j]]<-KLsALL
KLsALL<-rbind(KLsALL,sum(KLs))
KLsALLRandom<-rbind(KLsALLRandom,sum(KLsRandom))
explainKL<- which(s)
results<-list(explainKL=explainKL,KLsALL=KLsALL, KLsALLRandom)
return(results)
}


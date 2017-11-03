KLdivSelectVars<-function(d1,d2,m){
  eigenvalues1 <- (d1$sdev)^2
  eigenvalues1<-eigenvalues1[m]
  eigenvectors1 <-d1$rotation[m,]
  eigenvalues2 <- (d2$sdev)^2
  eigenvalues2<-eigenvalues2[m]
  eigenvectors2 <-d2$rotation[m,]
  
  cP<-min(which(stats:::summary.prcomp(d1)$importance[3, ]>=0.5))
  cQ<-min(which(stats:::summary.prcomp(d2)$importance[3, ]>=0.5))
  a<-sum(eigenvalues1[1:(cP-1)])
  b<-0.5*nrow(eigenvectors1)
  c<-sum(eigenvalues1[1:cP])
  if ( (c>=b) & (b>a) ){
    check<-TRUE 
  } else {
    check<-FALSE 
  }
  
  d<-sum(eigenvalues2[1:(cQ-1)])
  e<-0.5*nrow(eigenvectors2)
  f<-sum(eigenvalues2[1:cQ])
  if ( (f>=e) & (e>d) ){
    check2<-TRUE 
  } else {
    check2<-FALSE 
  }
  
  
  eigenvalues1R<-c( ((b/c)* eigenvalues1[1:cP] + 0.5) ,  rep(0.5, (nrow(eigenvectors1)-cP)) )
  eigenvalues2R<-c( ((e/f)* eigenvalues2[1:cQ] + 0.5) ,  rep(0.5, (nrow(eigenvectors2)-cQ)) )
  
  s2<-c()
  for (i in 1:cP){
    s1<-c()
    for (j in 1:cQ){
      s1[j]<-((eigenvalues1R[i]/0.5)-1)*( 1- (0.5/eigenvalues2R[j])) * 
        (eigenvectors1[,i]%*%eigenvectors2[,j])^2
    }
    s2[i]<-sum(s1)
  }
  
  A<-sum(s2)
  
  s2<-c()
  for (i in 1:cP){
    s1<-c()
    for (j in 1:cQ){
      s1[j]<-( 1- (0.5/eigenvalues1R[i])) *((eigenvalues2R[j]/0.5)-1)* 
        (eigenvectors1[,i]%*%eigenvectors2[,j])^2
    }
    s2[i]<-sum(s1)
  }
  B<-sum(s2)
  
  C<-sum(eigenvalues1R[1:cP]/0.5)+sum(eigenvalues2R[1:cQ]/0.5)
  D<-sum(0.5/eigenvalues1R[1:cP])+sum(0.5/eigenvalues2R[1:cQ])
#  MinSKL<-sum( (eigenvalues1R/eigenvalues2R) + (eigenvalues2R/eigenvalues1R) )- 2*nrow(eigenvectors1)   
 # MaxSKL<-sum( (eigenvalues1R/rev(eigenvalues2R)) + (rev(eigenvalues2R)/eigenvalues1R) )- 2*nrow(eigenvectors1)   
  
  
  SKL<-0.5*(C+D-(2*(cP+cQ))-A-B)
  
  return(SKL)
}
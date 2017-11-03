compareDsets<-function(d1,d2){
  a<-mat.mat(d1,d2)
  b<-which(a==1)
  return(length(b))
}

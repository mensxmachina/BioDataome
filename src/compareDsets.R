#this function receives as input two datasets d1 and d2
# and returns how many samples are shared between those datasets
#the comparison is column-wise


compareDsets<-function(d1,d2){
  a<-mat.mat(d1,d2)
  b<-which(a==1)
  return(length(b))
}

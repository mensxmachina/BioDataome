#find the average KL of a clique
#dsets is a character vector of dsets
#i.e dsets<-c("GSE26639","GSE54002","GSE20685")
#df is a numeric square symmetric numeric matrix of KLs 


cliqueKL<-function(dsets,df){
 kl<-c()
  combs<-t(combn(dsets,2))
   for (j in 1:nrow(combs)){
      a<-which(rownames(df)==combs[j,1])
      b<-which(colnames(df)==combs[j,2]) 
     kl[j]<-df[a,b]
   }
   kls<-mean(kl)
return(kls)
}
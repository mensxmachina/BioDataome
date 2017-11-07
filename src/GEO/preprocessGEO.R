#This function calls the SCAN method as described in 
#Piccolo SR, Sun Y, Campbell JD, Lenburg ME, Bild AH and Johnson WE (2012). 
#“A single-sample microarray normalization method to facilitate personalized-medicine workflows.” 
#Genomics, 100(6), pp. 337-344.
#it receives as input the path where the CEL files are stored and the number of cores to 
#run in parallel and returns a matrix of dimension: probes x samples with the normalized expression values

preprocessGEO<-function(x,y){
  setwd(x)
  cl <- makeCluster(y)
  registerDoParallel(cl)
  normalized =SCAN( "*.CEL.gz")
  stopCluster(cl)
  normalized<-exprs(normalized)
  return(normalized)
  
}



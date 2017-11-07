#this function receives as input a GEO Series (GSE) id, i.e."GSE11761"
#and returns a list with all series matrices. The same GSE may be linked 
#to more than one platforms (i.e GPL570 and GPL1261). 
#The size of the output list is the number of the measuring platforms. 

downloadPhenotype <- function(x) {
  phenotype <- NULL
  a <- TRUE
  while(a) {
    phenotype <- try(getGEO(x, GSEMatrix=TRUE))
    if(!inherits(phenotype,"try-error")) a <- FALSE
    else cat("NCBI was busy, trying again\n")
  }
  return(phenotype)
}

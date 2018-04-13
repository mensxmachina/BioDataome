#' Annotate a study (GSE) with a disease term from the Disease Ontology by exploiting only GEO
#'
#' Given a GSE id this function annotates the study with a disease term from the Disease Ontology (D-O):
#' http://disease-ontology.org/
#' It provides the most specific disease term, meaning the term with the highest depth in the D-O.
#'
#' @param GSE a GSE series ID
#' @return a character vector of all related diseases, separated by ;
#' @examples
#' diseases<-GSEtoDiseaseGEO("GSE10245")
#' @export

GSEtoDiseaseGEO<-function(GSE){
  if (missing(GSE))
    stop("Need to specify a GEO Series id, i.e 'GSE10026'")
  if (!grepl("GSE[0-9]+",GSE))
    stop("GSE must be a GEO Series id, i.e 'GSE10026'")

  #unlist all gses per disease
  tt<-strsplit(diseasesLevel$GSEs,";")
  #isolate numerical part of GSE code
  gse<-strsplit(GSE,"GSE")
  gse<-unlist(gse)[[2]]
  #look for GSE in gses per disease list
  c<-mapply(match,gse,tt)
  d<-which(!(is.na(c)))
  #if this GSE is found
  if (length(d)>0){
    diseases<-diseasesLevel[d,1]
    #remove "disease" which causes NAs
    diseases<-diseases[which(diseases!="disease")]
    if (length(diseases)>0){
      #find all leafs
      leaf<-c()
      for (j in 1:length(diseases)){
        leaf[j]<-match(tolower(diseases[j]),tolower(leafs[,1]))
      }
      #if we have only one leaf keep that
      if (length(which(!is.na(leaf)))==1){
        diseases<-diseases[which(!is.na(leaf))]
      } else if (length(which(!is.na(leaf)))>1){
        #if we have more than one leaves find parent categories of all entities
        #and keep the leaf that belongs to the most common category
        aa<-diseaseSubCategoryALLU[match(tolower(diseases),tolower(diseaseSubCategoryALLU[,1])),4]
        commonCategory<-names(table(aa)[which(table(aa)==max(table(aa)))])
        leafCats<-aa[which(!is.na(leaf))]
        #which leafs belong to the most common category
        leafsToKeep<-which(leafCats %in% commonCategory)
        #if no leaf belongs to the common category keep them all
        if (length(leafsToKeep)==0){
          diseases<-paste(diseases[which(!is.na(leaf))], sep=";", collapse=";")
        } else {
          diseases<-paste(diseases[leafsToKeep], sep=";", collapse=";")
        }
      } else {
        #this means that all diseases are nodes
        #find deepest nodes
        depth<-c()
        for (j in 1:length(diseases)){
          depth[j]<-as.numeric(diseaseSubCategoryALLU[match(tolower(diseases[j]),tolower(diseaseSubCategoryALLU[,1])),3 ] )
        }
        #if we have more than one category with the same depth, keep the most representative as in leafs

        bb<-which(depth==max(depth,na.rm = T))
        #if we have only one deepest node keep that
        if (length(bb)==1){
          diseases<-diseases[which(depth==max(depth,na.rm = T))]
        } else if (length(bb)>1){
          aa<-diseaseSubCategoryALLU[match(tolower(diseases),tolower(diseaseSubCategoryALLU[,1])),4]
          commonCategory<-names(table(aa)[which(table(aa)==max(table(aa)))])
          nodesToKeep<-which(aa[bb] %in% commonCategory)
          if (length(nodesToKeep)==0){
            diseases<-paste(diseases[bb], sep=";", collapse=";")
          } else {
            diseases<-paste(diseases[bb[nodesToKeep]], sep=";", collapse=";")
          }
        }
      }
    } else {
      diseases<-NA
    }
  } else {
    diseases<-NA
  }
  return(diseases)
}

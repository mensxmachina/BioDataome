#' Map a Disease Ontology (D-O) term to the first children nodes in D-O
#'
#' This function uses internal look up data to map a disease to its first children node.
#'
#' @param x a disease in D-O terms
#' @return the first children node of x disease
#' @examples
#' DOChild<-diseasetoChildrenNodes("vesiculitis")
#' @export

diseasetoChildrenNodes<-function(x){
  DOChild<-BioDataome:::diseaseSubCategoryALLU[match(tolower(x),tolower(BioDataome:::diseaseSubCategoryALLU[,1])),4]
  return(DOChild)
}

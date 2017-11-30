#' Map a Disease Ontology (D-O) term to the parent nodes in D-O
#'
#' This function uses internal look up data to map a disease to its parent node.
#'
#' @param x a disease in D-O terms
#' @return the parent node of x disease
#' @examples
#' DOParent<-diseasetoParentNodes("vesiculitis")
#' @export

diseasetoParentNodes<-function(x){
  DOParent<-BioDataome:::diseaseCategoryALLU[match(tolower(x),tolower(BioDataome:::diseaseCategoryALLU[,1])),4]
  return(DOParent)
}


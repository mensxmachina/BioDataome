#' Discover control samples from phenotype data in GEO
#'
#' This function discovers control samples from the series matrix found in GEO
#' It searches for specific keywords that are often used to denote controls, in specific columns of
#' series matrices.
#'
#' @param x a data frame with the contents of series matrix
#' @return a data frame of GEO sample ids (i.e. GSM60555) and their class.
#' @examples
#' phenos<-downloadPhenotypePlatform("GSE11761","GPL570")
#' controls<-controlSamples(phenos)
#' @export


controlSamples<-function(d){

  controls <- c("untreated", "control", "ctrl","healthy",
               "unstimulated", "unaffected", "normal", "non","undosed")
  TitleControls <- grep(paste(controls,collapse="|"),
                       as.character(d$title), ignore.case = T)
  SourceControls <- grep(paste(controls,collapse="|"),
                        as.character(d$source_name_ch1), ignore.case = T)
  CharacterControls <- grep(paste(controls,collapse="|"),
                           as.character(d$characteristics_ch1), ignore.case = T)
  #merge control samples found in any of the three columns of series matrix (title, source_name_ch1, or characteristics_ch1 )
  controlsU<-c(rownames(d)[TitleControls],rownames(d)[SourceControls],rownames(d)[CharacterControls])
  #create a data frame of sample names and class
  class<-character(length = nrow(d))
  class[match(unique(controlsU),rownames(d))]<-"control"
  class[setdiff(1:nrow(d),match(unique(controlsU),rownames(d)))]<-"unknown"

  controls<-cbind(GSM=rownames(d),class=class)

  return(controls)

}

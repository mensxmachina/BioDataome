#' Run all steps to download, preprocess and annotate a GEO dataset for BioDataome
#'
#' Given a GSE id this function downloads, preprocesses, annotates a study and also creates the metadata. It saves two files,
#' the data file and the GEO metadata file in the given path and returns a character vector of BioDataome metadata with
#' with information for GSE, Species,	Entity,	Technology,	Type,	number of samples, dublicate GSEs,	Disease,	DOParent,
#' DOChild, path to csv, path to D-O, path to GEO, path to metadata
#' @param x a GSE series ID
#' @param y a GEO platform id (GPL)
#' @param z the path to write the output
#' @return writes in the given path the preprocessed data, the metadata file and two data frames a character vector of metadata for GSE,	Species,	Entity,	Technology,	Type,	number of samples,
#' 	dublicates,	Disease,	DOParent,	DOChild, path to data, path to D-O, path to GEO, path to metadata

#' @examples
#' diseases<-GSEtoDisease("GSE10245")
#' @export
#' @importFrom rentrez entrez_search
#' @importFrom RCurl getURL
#' @importFrom XML xmlToList

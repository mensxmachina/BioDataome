#' Annotate a study (GSE) with a disease term from the Disease Ontology by exploiting both PubTator and GEO
#'
#' Given a GSE id this function annotates the study with a disease term from the Disease Ontology (D-O):
#' http://disease-ontology.org/
#' It provides the most specific disease term, meaning the term with the highest depth in the D-O.
#'
#' @param GSE a GSE series ID
#' @return a character vector of all related diseases, separated by ;
#' @examples
#' diseases<-GSEtoDisease("GSE10245")
#' @export
#' @importFrom rentrez entrez_search
#' @importFrom RCurl getURL
#' @importFrom XML xmlToList
#' @importFrom pubmed.mineR pubtator_function
#' @importFrom httr GET stop_for_status content

GSEtoDisease<-function(GSE){
  ### Intergrate Local PubTator

  if (missing(GSE))
    stop("Need to specify a GEO Series id, i.e 'GSE10026'")
  if (!grepl("GSE[0-9]+",GSE))
    stop("GSE must be a GEO Series id, i.e 'GSE10026'")

  #unlist all gses per disease
  tt<-strsplit(diseasesLevel$GSEs,";")

  #find Pubmed IDs for this GSE
  a<-paste0(GSE," [ACCN] AND gse[ETYP]")
  r_search<- rentrez::entrez_search(db="gds", term=a,retmax =10000, use_history=TRUE)

  if (r_search$count!=0){

    b<-entrez_summary(db="gds", id=r_search$ids[1])
    pubmed<-b$pubmedids

    if (length(pubmed)!=0){

      #use PubTator's API to annotate gse
      diseases<-c()

      for (j in 1:length(pubmed)){
        #uri<-paste("https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/Disease/",pubmed[j],"/BioC",sep="") ### Error
        #uri<-paste("https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/BioConcept/",pubmed[j],"/PubTator",sep="")
        uri<-paste0("https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids=",pubmed[j],"&concepts=disease")

        ###### ------ Error Tests ------ #####

        # opt <- RCurl::curlOptions(verbose = TRUE, sslversion = CURL_SSLVERSION_TLSv1_2)
        #
        # RCurl::getURL(uri, .opts = opt)
        #
        # uri <- "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids=28483577&concepts=disease"
        # uri <- "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml?pmcids=PMC6207735"
        #
        # require(RCurl)
        # #Works, note I removed sslversion
        # RCurl::getURL(uri, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"), .opts=c(followLocation=T, verbose = TRUE, sslversion=6L), ssl.verifypeer = FALSE)
        #
        # RCurl::getURL(uri, ssl.verifypeer = F,sslversion=6L)
        #
        # require(rvest)
        # t <- read_html(uri)
        #
        # require(curl)
        # t <- curl(uri)
        #
        # require(curl)
        # con <- curl(uri)
        # readLines(con)

        ###### ---- Error Test ---------- #######

        #require(httr)
        req <- GET(uri)
        stop_for_status(req)
        t <- content(req)
        diseases <- c(diseases,tolower(strsplit(t,"\t")[[1]][4]))

        ### Erroniously Block
        RCurl.attempt <- try(t <- RCurl::getURL(uri, ssl.verifypeer = F,sslversion=6L, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")), T)

        if (!(isTRUE(class(RCurl.attempt) == "try-error"))){

          if( !is.na(pmatch('[Error]', t)) ){

            diseases<-NA

          } else {

            ll <- XML::xmlToList(t)
            list1 <- ll$document[[3]]
            if (length(list1)>3){
              d <- c()
              for (k in 1:(length(list1)-3)){
                d[k] <- list1[[k+3]]$text
              }
              diseases <- c(diseases,tolower(d))
            }
          }
        }

        #####

        #------------
        #local.query <- getPubTatorQuery(pubmed[j])
        #local.query <- strsplit(local.query$MENTIONS,split='|', fixed=TRUE)
        #local.query <- lapply(local.query, function(x) tolower(x[which.max(nchar(x))]))
        #local.query <- unlist(local.query)
        #diseases <- c(diseases,local.query)
        #-----------

        pubtator.attemp <- try(local.query <- pubtator_function(pubmed[j]), T)

        if (!(isTRUE(class(pubtator.attemp)=='try-error'))){
          if (local.query !=" No Data ") {
            local.query <- local.query$Diseases
            if (!is.null(local.query)){
              local.query <- strsplit(local.query,split='>', fixed=TRUE)
              local.query <- lapply(local.query, function(x) tolower(x[1]))
              diseases <- c(diseases, unlist(local.query))
            }

          }
        }
        #-----------
        }

      }

      diseases<-unique(diseases)
      #if PubTator returns an empty set
      #try with GEO
      if (length(diseases)==0){

        diseases<-sapply(GSE,GSEtoDiseaseGEO)

      } else {

        #find all leafs
        leaf<-c()

        for (j in 1:length(diseases)){

          leaf[j]<-match(tolower(diseases[j]),tolower(leafs[,1]))

        }

        #if we have only one leaf keep that
        if (length(which(!is.na(leaf)))==1){

          diseases<-diseases[which(!is.na(leaf))]

        } else if (length(which(!is.na(leaf)))>1){

          #if we have more than one leafs find parent categories of all entiites
          #and keep the leaf that belongs to the most common category

          aa <- diseaseSubCategoryALLU[match(tolower(diseases),tolower(diseaseSubCategoryALLU[,1])),4]

          commonCategory <- names(table(aa)[which(table(aa)==max(table(aa)))])

          leafCats <- aa[which(!is.na(leaf))]

          #which leafs belong to the most common category
          leafsToKeep <- which(leafCats %in% commonCategory)

          #if no leaf belong to the common category keep them all

          if (length(leafsToKeep)==0){

            diseases<-paste(diseases[which(!is.na(leaf))], sep=";", collapse=";")

          } else {

            diseases<-paste(diseases[match(leafCats[leafsToKeep],aa)], sep=";", collapse=";")

          }
        } else {

          #if we only have nodes
          #find deepest nodes
          depth<-c()
          for (j in 1:length(diseases)){

            depth[j]<-as.numeric(diseaseSubCategoryALLU[match(tolower(diseases[j]),tolower(diseaseSubCategoryALLU[,1])),3 ] )####

          }
          #if we have more than one category with the same depth, keep the most representative as in leafs

          bb<-which(depth==max(depth,na.rm = T))

          if (length(bb)==0){

            diseases<-sapply(GSE,GSEtoDiseaseGEO)

          } else if (length(bb)==1){

            diseases<-diseases[which(depth==max(depth,na.rm = T))]

          } else {

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
      }

      if (is.na(diseases)){

        diseases<-sapply(GSE,GSEtoDiseaseGEO)

      }


    } else {
      diseases<-sapply(GSE,GSEtoDiseaseGEO)
    }

  # } else {
  #   diseases<-NA
  # }

  return(diseases)
}



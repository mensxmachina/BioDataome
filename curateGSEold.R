#' Run all steps to download, preprocess and annotate a GEO dataset
#'
#' Given a GSE id this function downloads, preprocesses, annotates a study and also creates the sample phenotype metadata.

#' @param x a GSE series ID
#' @param y a GEO technology id (GPL)
#' @param z the path to save the downloaded files. By default this value is set to the working directory.
#' @return a list with three components:
#' \itemize{
#'  \item{"metadata"}{a data frame of the metadata with sample phenotype information}
#'  \item{"dataNorm"}{a data frame of the preprocessed data with variables (probes) as rows and samples as columns}
#'  \item{"indexdata"}{a character vector with all the necessary information that describe the dataset}
#' }
#' @examples
#' curateGSE("GSE11761","GPL570",getwd())
#' @importFrom GEOquery gunzip
#' @export

  curateGSE<-function(x,y,z=getwd(), keepRaw=FALSE){
    #platforms currenty curated in BioDataome
    platforms<-c("GPL570", "GPL96", "GPL6244", "GPL1261","GPL13534")

    if (missing(x))
      stop("Need to specify a GEO Series id, i.e 'GSE10026'")
    if ( (missing(y)) | (!file.exists(y)) )
      stop("Need to specify a GEO technology id (GPL570)")
    if (!grepl("GSE[0-9]+",x))
      stop("x must be a GEO Series id, i.e 'GSE10026'")
    if (!any(grepl(y, platforms, ignore.case=TRUE)))
      stop("y must be one of the technologies: 'GPL570', 'GPL96', 'GPL6244', 'GPL1261','GPL13534'")

    #download and curate phenotype metadata
    metadata<-GSEmetadata(x,y)
    ##save(metadata, file=paste0(z,"/",x,"_Annot.Rda"))
    #download raw files
    downloadRaw(x,z)
    untar(file.path(z,x,paste0(x,"_RAW.tar")), exdir = file.path(z,x))
    #remove compressed file
    file.remove(file.path(z,x,paste0(x,"_RAW.tar")))
    #by default downloadRaw will download all files for the specific GSE

    #if the user specified GPL13534 with is the Illumina HumanMethylation450 BeadChip
    #list and preprocess idat files. Else list and preprocess CEL files

    if (y=="GPL13534"){
      idatFiles <- list.files(file.path(z,x), pattern = "idat.gz$", full = TRUE)
      sapply(idatFiles, gunzip, overwrite = TRUE)
      idatFiles2 <- list.files(file.path(z,x), pattern = "idat$",include.dirs = FALSE)
      #keep only samples found in phenos
      idatKeep<-c()
      for (j in 1:nrow(metadata)){
        if (any(startsWith(idatFiles2,as.character(metadata$samples[j])))==TRUE){
          idatKeep<-c(idatKeep,which(startsWith(idatFiles2,as.character(metadata$samples[j]))))
        }
      }
      idatFiles2<-idatFiles2[idatKeep]
      #replace lower grn with upper
      replace_grn<-function(x){
        return(gsub(pattern = "grn", replacement = "Grn", x))
      }
      #replace lower red with upper
      replace_red<-function(x){
        return(gsub(pattern = "red", replacement = "Red", x))
      }

      sall<-unlist(lapply(idatFiles2, replace_grn))
      sall<-unlist(lapply(sall, replace_red))

      basenames<-unlist(strsplit(sall, "_Grn.idat"))
      basenames<-unique(unlist(strsplit(basenames, "_Red.idat")))

      basenames<-file.path(file.path(z,x), basenames)
      dataNorm<-preprocessGEOMethylation(x)

    } else {

      celfiles<-list.files(file.path(z,x), pattern = ".cel",ignore.case = T)
      #keep only cel files that have annotation
      celsKeep<-celfiles[startsWith(celfiles,as.character(metadata$samples))]
      #remove the rest
      toRemove<-setdiff(celfiles,celsKeep)
      file.remove(file.path(z,x,toRemove))
      dataNorm<-preprocessGEO(file.path(z,x),3)

    }

    #save Rda file
    save(dataNorm, file=paste0(z,"/",x,".Rda"))
    #remove directory with raw files unless the user specifies keepRaw=TRUE
    if (keepRaw==FALSE){
      file.remove(dir(file.path(z,x), full.names=TRUE))
    }
  }
